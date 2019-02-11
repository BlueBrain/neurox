#include "neurox/interpolators/variable_time_step.h"

#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver */
#include "cvode/cvode.h"

#include "cvode/cvode_spils.h"

#include <set>

using namespace neurox;
using namespace neurox::interpolators;
using namespace neurox::tools;

const char *VariableTimeStep::GetString() { return "VariableTimeStep"; }

void VariableTimeStep::CopyState(Branch *branch, N_Vector y, const CopyOp op) {
  const VariableTimeStep *vardt = (VariableTimeStep *)branch->interpolator_;
  const bool is_scatter_op =
      op == CopyOps::kScatterYdot || op == CopyOps::kScatterY;
  const bool use_ydot =
      op == CopyOps::kScatterYdot || op == CopyOps::kGatherYdot;
  double *y_data = NV_DATA_S(y);  // y or ydot
  double **var_map = use_ydot ? vardt->state_dv_map_ : vardt->state_var_map_;

  const int &cap_count =
      branch->mechs_instances_[mechanisms_map_[CAP]].nodecount;
  const int iters_limit = LAYOUT == 0 ? cap_count : vardt->equations_count_;

  // AoS will map all variables, SoA only map no_cap voltages/RHS for now
  int i = -1;
  if (is_scatter_op)
    for (i = 0; i < iters_limit; i++) *(var_map[i]) = y_data[i];
  else
    for (i = 0; i < iters_limit; i++) y_data[i] = *(var_map[i]);

#if LAYOUT == 0
  // SoA mapping takes advantage of state vars sequential mem-alignment
  for (int m = 0; m < neurox::mechanisms_count_; m++) {
    const int state_vars_count = neurox::mechanisms_[m]->state_vars_->count_;
    const int nodecount = branch->mechs_instances_[m].nodecount;
    if (is_scatter_op)
      for (int s = 0; s < state_vars_count; s++, i += nodecount)
        memcpy(var_map[i], &y_data[i], sizeof(double) * nodecount);
    else
      for (int s = 0; s < state_vars_count; s++, i += nodecount)
        memcpy(&y_data[i], var_map[i], sizeof(double) * nodecount);
  }
#endif
}

void VariableTimeStep::ScatterY(Branch *branch, N_Vector y) {
  CopyState(branch, y, CopyOps::kScatterY);
}

void VariableTimeStep::GatherY(Branch *branch, N_Vector y) {
  CopyState(branch, y, CopyOps::kGatherY);
}

void VariableTimeStep::ScatterYdot(Branch *branch, N_Vector ydot) {
  CopyState(branch, ydot, CopyOps::kScatterYdot);
}

void VariableTimeStep::GatherYdot(Branch *branch, N_Vector ydot) {
  CopyState(branch, ydot, CopyOps::kGatherYdot);
}

int VariableTimeStep::RootFunction(realtype t, N_Vector y, realtype *gout,
                                   void *user_data) {
  Branch *branch = (Branch *)user_data;

  // get offset and voltage of Axon Initial Segment
  int ais_offset = branch->thvar_ptr_ - branch->nt_->_actual_v;
  double v_ais = NV_Ith_S(y, ais_offset);

  // How it works: when gout[x] is zero, a root is found
  // assert(v_ais >= -150 && v_ais < 50);
  gout[0] = v_ais - branch->soma_->threshold_;  // AP threshold reached
  return CV_SUCCESS;
}

int VariableTimeStep::RHSFunction(realtype t, N_Vector y, N_Vector ydot,
                                  void *user_data) {
  Branch *branch = (Branch *)user_data;
  VariableTimeStep *vardt = (VariableTimeStep *)branch->interpolator_;
  NrnThread *nt = branch->nt_;

  //////// occvode.cpp: Cvode::fun_thread_transfer_part1

  printf("\nBRUNO BEFORE t=%f fun_thread\n", nt->_t);
  for (int i = 0; i < NV_CONTENT_S(y)->length; i++)
    printf("BRUNO y[%d]=%f\n", i, NV_CONTENT_S(y)->data[i]);
  if (ydot)
    for (int i = 0; i < NV_CONTENT_S(ydot)->length; i++)
      printf("BRUNO y'[%d]=%f\n", i, NV_CONTENT_S(ydot)->data[i]);

  const double h = vardt->cvode_mem_->cv_h;
  nt->_dt = h == 0 ? 1e-8 : h;  // TODO set 1e-8 to input_params_->dt/2
  nt->cj = 1 / nt->_dt;
  nt->_t = t;

  // Updates internal states of continuous point processes (vecplay)
  // e.g. stimulus. vecplay->pd points to a read-only var used by
  // point proc mechanisms' nrn_current function
  // cvtrset.cpp :: CVode::fun_thread_transfer_part1()

  branch->FixedPlayContinuous(nt->_t);

  // copies V and state-vars from CVODES to NrnThread
  VariableTimeStep::ScatterY(branch, y);

  // start of occvode.cpp :: nocap_v
  HinesSolver::ResetRHSandDNoCapacitors(branch);

  // sum mech-instance contributions to D and RHS on no-caps
  branch->CallModFunction(Mechanism::ModFunctions::kCurrent,
                          vardt->no_cap_ml_);  // rhs
  branch->CallModFunction(Mechanism::ModFunctions::kJacob,
                          vardt->no_cap_ml_);  // lhs

  HinesSolver::SetupMatrixVoltageNoCapacitors(branch);

  //////// ocvode2.cpp: Cvode::fun_thread_transfer_part2

  // cvtrset.cpp :: CVode::rhs
  HinesSolver::ResetArray(branch, nt->_actual_rhs);

  // cvtrset.cpp :: CVode::rhs() -> rhs_memb(z.cv_memb_list_)
  // sum mech-instance contributions to D and RHS on caps
  branch->CallModFunction(Mechanism::ModFunctions::kCurrent);

  // add parent and children axial currents (A*dv and B*dv) to RHS
  // cvtrset.cpp :: CVode::rhs()
  HinesSolver::SetupMatrixRHS(branch);

  // update mechanisms state (eg opening vars and derivatives)
  // cvtrset.cpp :: CVode::fun_thread_transfer_part2() -> do_ode()
  branch->CallModFunction(Mechanism::ModFunctions::kODESpec);

  // divide RHS by Cm and compute capacity current
  // cvtrset.cpp :: CVode::fun_thread_transfer_part2() -> nrn_div_capacity()
  branch->CallModFunction(Mechanism::ModFunctions::kDivCapacity);

  // copies dV/dt (RHS) and state-vars-derivative to CVODES
  VariableTimeStep::GatherYdot(branch, ydot);

  printf("\nBRUNO AFTER t=%f fun_thread\n", nt->_t);
  for (int i = 0; i < NV_CONTENT_S(y)->length; i++)
    printf("BRUNO y[%d]=%f\n", i, NV_CONTENT_S(y)->data[i]);
  if (ydot)
    for (int i = 0; i < NV_CONTENT_S(ydot)->length; i++)
      printf("BRUNO y'[%d]=%f\n", i, NV_CONTENT_S(ydot)->data[i]);

  return CV_SUCCESS;
}

/* The solve function must solve the linear system Mx=b, where M is an
approximation of the Newton matrix, I - gamma J, where J = df/dy,
and the rhs-vector b is an input. Called once per newton, thus
several times per time-step. (occvode.cpp: Cvode::solvex_thread) */
int VariableTimeStep::PreConditionedDiagonalSolver(CVodeMem cv_mem, N_Vector b,
                                                   N_Vector weight,
                                                   N_Vector ycur,
                                                   N_Vector fcur) {
  // b is the right-hand-side vector, solution to be returned in b
  // ycur contains vector approximations to y(t_n)
  // ycur contains vector approximations to f(t_n, ycur)
  Branch *branch = (Branch *)cv_mem->cv_user_data;
  NrnThread *nt = branch->nt_;
  nt->_dt = cv_mem->cv_gamma;
  nt->cj = 1.0 / nt->_dt;

  // Cvode::lhs()
  HinesSolver::ResetArray(branch, nt->_actual_d);
  branch->CallModFunction(Mechanism::ModFunctions::kJacob);
  branch->CallModFunction(Mechanism::ModFunctions::kJacobCapacitance);
  assert(nt->_actual_b[0] == 0);
  HinesSolver::SetupMatrixDiagonal(branch);
  // end of Cvode::lhs()

  ScatterYdot(branch, b);
  branch->CallModFunction(Mechanism::ModFunctions::kMulCapacity);
  HinesSolver::ResetRHSNoCapacitors(branch);
  HinesSolver::BackwardTriangulation(branch);
  HinesSolver::ForwardSubstitution(branch);
  branch->CallModFunction(Mechanism::ModFunctions::kODEMatsol);
  GatherYdot(branch, b);
  return CV_SUCCESS;
}

// jacobian routine: compute J(t,y) = df/dy
int VariableTimeStep::JacobianDense(realtype t, N_Vector y, N_Vector fy,
                                    SUNMatrix J, void *user_data, N_Vector,
                                    N_Vector, N_Vector) {
  _SUNMatrixContent_Dense *content = (_SUNMatrixContent_Dense *)J->content;
  realtype **jac = content->cols;
  Branch *branch = (Branch *)user_data;
  VariableTimeStep *vardt = (VariableTimeStep *)branch->interpolator_;
  NrnThread *nt = branch->nt_;
  assert(t == nt->_t);

  // RHS provided the righ-hand side
  // We will now compute the LHS
  // Newton Matrix: M = I - gamma*J

  // in neuron we solve Px=b with P approximates Id-gamma*J
  // const double gamma = vardt->cvode_mem_->cv_gamma;
  // nt->_dt = gamma;
  // nt->cj = 1/nt->_dt;

  // used in kJacobCapacitance
  nt->_dt = 1;
  nt->cj = 1;

  // if I reset V, then nrn_current is wrong!
  // cvtrset.cpp :: CVode:: lhs
  // HinesSolver::ResetMatrixV(branch);

  // does nothing so far (called before as part of nrn_current)
  // cvtrset.cpp :: CVode:: lhs -> lhs_memb()
  branch->CallModFunction(Mechanism::ModFunctions::kJacob);

  // D = D + cfac * .001 * _nt->cj // decay of D?
  branch->CallModFunction(Mechanism::ModFunctions::kJacobCapacitance);

  /*
  //add axial currents to D (d[i] -= b[i] and d[p[i]] -= a[i])
  HinesSolver::SetupMatrixDiagonal(branch);

  //RHS = RHS* (.001 * nt->cj;)*cm; //decay of dV/dt
  branch->CallModFunction(Mechanism::ModFunctions::kMulCapacity);

  branch->SolveTreeMatrix(); //Gaussian Elimination
  //now we have 0A + 1D + 0B = RHS, ie dV/dt = RHS

  //if mechs stiff()==2.
  //branch->CallModFunction(Mechanism::ModFunctions::kODEMatsol);
  */
  /// end of neuron occvode.cpp::solvex_thread()

  const int capacitors_count = nt->end;
  const double *a = nt->_actual_a;
  const double *b = nt->_actual_b;
  const double *d = nt->_actual_d;
  const int *p = nt->_v_parent_index;

  // Jacobian for main current equation:
  // d(dV/dt) / dV     = sum_i g_i x_i  + D  //mechs currents + D
  // d(dV/dt) / dV_p   = -A          //parent compartment
  // d(dV/dt) / dV_c_i = -B_c_i  //children_i compartment

  for (int n = 0; n < capacitors_count; n++) {
    jac[n][n] = d[n];  // D = d (dV_n/dt) /dV_n
    if (n == 0) continue;
    jac[p[n]][n] = a[n];  // A = d (dV_p/dt) /dV_n
    jac[n][p[n]] = b[n];  // B = d (dV_n/dt) /dV_p
  }

  // Matsol  Matsol solves (1 + dt*jacobian)*x = b
  // for direct diagonal matrix solver
  // branch->CallModFunction(Mechanism::ModFunctions::kODEMatsol);

  // int compartment_id=-1;
  int dv_offset = capacitors_count;
  for (int m = 0; m < neurox::mechanisms_count_; m++) {
    Mechanism *mech = mechanisms_[m];
    Memb_list *mech_instances = &branch->mechs_instances_[m];
    for (int n = 0; n < mech_instances->nodecount; n++)
      for (int s = 0; s < mech->state_vars_->count_; s++) {
        // Reminder: state vars Dm represents d(dm_i/dt) / dm_j
        jac[dv_offset][dv_offset] =
            *(vardt->state_dv_map_[dv_offset - capacitors_count]);
        dv_offset++;
      }
  }
  printf("Jac neuron %d, t%.10f\n", branch->soma_->gid_, t);
  double *y_data = NV_DATA_S(y);
  // for (int i=0; i<N; i++)
  //    printf("t=%f\t RHS[%d] = %.10f\t Jac[%d] = %.10f\n", t, i, y_data[i], i,
  //    jac[i][i]);

  return CV_SUCCESS;
}

VariableTimeStep::VariableTimeStep()
    : cvode_mem_(nullptr),
      equations_count_(-1),
      state_var_map_(nullptr),
      state_dv_map_(nullptr),
      y_(nullptr),
      no_cap_child_ids_(nullptr),
      no_cap_node_ids_(nullptr),
      no_cap_ml_(nullptr) {}

VariableTimeStep::~VariableTimeStep() {
  N_VDestroy_Serial(y_);
  CVodeFree((void **)(&cvode_mem_));

  delete no_cap_child_ids_;
  delete no_cap_node_ids_;
  // because these point to the same memory as
  // nt->data, we only delete main array
  delete[] no_cap_ml_;
}

// Neuron :: occvode.cpp :: init_global()
void VariableTimeStep::Init(Branch *branch) {
  assert(branch->interpolator_ != nullptr);
  VariableTimeStep *vardt = (VariableTimeStep *)branch->interpolator_;
  CVodeMem &cvode_mem = vardt->cvode_mem_;
  int cap_count = branch->mechs_instances_[mechanisms_map_[CAP]].nodecount;
  NrnThread *&nt = branch->nt_;

  int flag = CV_ERR_FAILURE;

  // some methods from Branch::Finitialize
  branch->CallModFunction(Mechanism::ModFunctions::kThreadTableCheck);
  branch->InitVecPlayContinous();
  branch->DeliverEvents(t);
  for (int n = 0; n < branch->nt_->end; n++)
    branch->nt_->_actual_v[n] = input_params_->voltage_;
  branch->CallModFunction(Mechanism::ModFunctions::kBeforeInitialize);
  branch->CallModFunction(Mechanism::ModFunctions::kInitialize);
  branch->CallModFunction(Mechanism::ModFunctions::kAfterInitialize);
  branch->CallModFunction(Mechanism::ModFunctions::kBeforeStep);
  branch->DeliverEvents(t);

  // equations: capacitors + mechanisms * states
  int &equations_count = vardt->equations_count_;
  equations_count = cap_count;
  for (int m = 0; m < neurox::mechanisms_count_; m++) {
    Mechanism *mech = mechanisms_[m];
    Memb_list *mech_instances = &branch->mechs_instances_[m];
    equations_count += mech_instances->nodecount * mech->state_vars_->count_;
  }

  ///// create initial state y_0 for state array y

  // create map from y and dy to NrnThread->data
  vardt->state_var_map_ = new double *[equations_count]();
  vardt->state_dv_map_ = new double *[equations_count]();

  int var_offset = 0;
  Memb_list *capac_instances = &branch->mechs_instances_[mechanisms_map_[CAP]];
  for (int c = 0; c < capac_instances->nodecount; c++) {
    int compartment_id = capac_instances->nodeindices[c];
    vardt->state_var_map_[var_offset] = &(nt->_actual_v[compartment_id]);
    vardt->state_dv_map_[var_offset] = &(nt->_actual_rhs[compartment_id]);
    var_offset++;
  }

  ////////// collect information about non-capacitors nodes

  vardt->no_cap_node_ids_count_ = nt->end - capac_instances->nodecount;
  vardt->no_cap_node_ids_ = new int[vardt->no_cap_node_ids_count_];
  int no_cap_count = 0;

  std::vector<int> child_ids;

  // get list of all nodes that are capacitors
  std::set<int> capacitor_ids;
  for (int c = 0; c < capac_instances->nodecount; c++) {
    int compartment_id = capac_instances->nodeindices[c];
    capacitor_ids.insert(compartment_id);
  }

  // inspired by neuron data structures
  for (int i = 0; i < nt->end; i++) {
    // if this node is not a capacitors node
    if (capacitor_ids.find(i) == capacitor_ids.end())
      vardt->no_cap_node_ids_[no_cap_count++] = i;

    // if parent node is not a capacitors node
    if (i > 0 &&
        capacitor_ids.find(nt->_v_parent_index[i]) == capacitor_ids.end())
      child_ids.push_back(i);
  }

  // create childs ids and count of no-cap parents
  vardt->no_cap_child_ids_count_ = child_ids.size();
  vardt->no_cap_child_ids_ = new int[vardt->no_cap_child_ids_count_];
  memcpy(vardt->no_cap_child_ids_, child_ids.data(),
         child_ids.size() * sizeof(int));
  assert(vardt->no_cap_node_ids_count_ == no_cap_count);

  // occvode.cpp::new_no_cap_memb()
  Vectorizer::GroupBranchInstancesByCapacitors(branch, &(vardt->no_cap_ml_),
                                               nullptr, &capacitor_ids);

  ////////////  build remaining map with state vars
  for (int m = 0; m < neurox::mechanisms_count_; m++) {
    Mechanism *mech = mechanisms_[m];
    Memb_list *mech_instances = &branch->mechs_instances_[m];
    int ml_data_offset = 0;

#if LAYOUT == 1
    for (int n = 0; n < mech_instances->nodecount; n++) {
      for (int s = 0; s < mech->state_vars_->count_; s++) {
#else
    for (int s = 0; s < mech->state_vars_->count_; s++) {
      for (int n = 0; n < mech_instances->nodecount; n++) {
#endif
        int state_var_index = mech->state_vars_->var_offsets_[s];
        int state_dv_index = mech->state_vars_->dv_offsets_[s];
#if LAYOUT == 1
        int state_var_offset =
            ml_data_offset + mech->data_size_ * n + state_var_index;
        int state_dv_offset =
            ml_data_offset + mech->data_size_ * n + state_dv_index;
#else
        // Note: We add SoA format to the maps, without padding because
        // it would require larger+padded y and y' arrays for CVODE
        int state_var_offset =
            ml_data_offset +
            Vectorizer::SizeOf(mech_instances->nodecount) * state_var_index + n;
        int state_dv_offset =
            ml_data_offset +
            Vectorizer::SizeOf(mech_instances->nodecount) * state_dv_index + n;
#endif
        assert(state_var_offset <
               Vectorizer::SizeOf(mech_instances->nodecount) *
                   mech->data_size_);
        assert(state_dv_offset < Vectorizer::SizeOf(mech_instances->nodecount) *
                                     mech->data_size_);
        vardt->state_var_map_[var_offset] =
            &(mech_instances->data[state_var_offset]);
        vardt->state_dv_map_[var_offset] =
            &(mech_instances->data[state_dv_offset]);
        var_offset++;
      }
    }
    ml_data_offset +=
        Vectorizer::SizeOf(mech_instances->nodecount) * mech->data_size_;
  }
  assert(var_offset == equations_count);
  vardt->y_ = N_VNew_Serial(equations_count);
  VariableTimeStep::GatherY(branch, vardt->y_);

  for (int i = 0; i < NV_CONTENT_S(vardt->y_)->length; i++)
    printf("BRUNO INIT y[%d]=%f\n", i, NV_CONTENT_S(vardt->y_)->data[i]);

  // absolute tolerance array (low for voltages, high for mech states)
  vardt->absolute_tolerance_ = N_VNew_Serial(equations_count);
  for (int i = 0; i < equations_count; i++) {
    double tol = i < cap_count ? input_params_->cvode_atol_v_
                               : input_params_->cvode_atol_states_;
    NV_Ith_S(vardt->absolute_tolerance_, i) = tol;
  }

  // Allocate mem block for BDF or Adams, with Newton solver (for stiff sol.)
  cvode_mem = (CVodeMem)CVodeCreate(CV_BDF, CV_NEWTON);

  // from cvodeobj.cpp :: cvode_init()
  vardt->cvode_mem_->cv_gamma = 0.;
  vardt->cvode_mem_->cv_h = 0.;

  // CVodeInit allocates and initializes memory for a problem
  double t0 = input_params_->tstart_;
  flag = CVodeInit(cvode_mem, VariableTimeStep::RHSFunction, t0, vardt->y_);
  assert(flag == CV_SUCCESS);

  // specify integration tolerances. MUST be called before CVode.
  flag = CVodeSVtolerances(cvode_mem, input_params_->cvode_rtol_,
                           vardt->absolute_tolerance_);
  assert(flag == CV_SUCCESS);

  // specify user data to be used on functions as void* user_data_ptr;
  flag = CVodeSetUserData(cvode_mem, branch);
  assert(flag == CV_SUCCESS);

  // specify root func. and roots (AP-threshold reached from below)
  int roots_direction[1] = {1};
  flag = CVodeRootInit(cvode_mem, 1, VariableTimeStep::RootFunction);
  CVodeSetRootDirection(cvode_mem, roots_direction);
  assert(flag == CV_SUCCESS);

  // Reminder: direct solvers give the solution (LU-decomposition, etc)
  // Indirect solvers require iterations (eg Jacobi method)
  switch (input_params_->interpolator_) {
    case InterpolatorIds::kCvodePreConditionedDiagSolver: {
      // CVODES guide chapter 8: Providing Alternate Linear Solver
      // Modules: only lsolve function is mandatory
      // (non-used functions need to be set to null)
      cvode_mem->cv_linit = nullptr;
      cvode_mem->cv_lsetup = nullptr;
      cvode_mem->cv_lfree = nullptr;
      // cvode_mem->cv_setupNonNull = FALSE;
      cvode_mem->cv_lsolve = PreConditionedDiagonalSolver;
      break;
    }
    case InterpolatorIds::kCvodeDenseMatrix: {
      /* Create dense SUNMatrix for use in linear solver */
      SUNMatrix A = SUNDenseMatrix(equations_count, equations_count);
      assert(A);

      /* Create dense SUNLinearSolver object for use by CVode */
      SUNLinearSolver LS = SUNDenseLinearSolver(this->y_, A);
      assert(LS);

      /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to
       * CVode */
      flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
      assert(flag);

      flag = CVDlsSetJacFn(cvode_mem, VariableTimeStep::JacobianDense);
      if (flag == CVDLS_MEM_FAIL) {
        throw std::runtime_error(
            "ERROR: can't allocate memory for dense jacobian for gid " +
            to_string(branch->soma_->gid_) + " and " +
            to_string(equations_count) + " equations.\n");
      }
      break;
    }
    case InterpolatorIds::kCvodeDiagonalMatrix:
      flag = CVDiag(cvode_mem);
      break;
    case InterpolatorIds::kCvodeSparseMatrix:
      assert(0);
      // Requires installation of Superlumt or KLU
      // flag = CVSlsSetSparseJacFn(cvode_mem, nullptr);
      // int nnz = equations_count * equations_count;
      // flag = CVKLU(cvode_mem, 1, equations_count, nnz);
      // flag = CVSuperLUMT(cvode_mem, 1, equations_count, nnz);

      // if not found, uses a dense matrix with sparse values (copy from above)
      break;
  }
  assert(flag == CV_SUCCESS);

  CVodeSetMinStep(cvode_mem, input_params_->dt_);
  CVodeSetMaxStep(cvode_mem, input_params_->tstop_);
  CVodeSetStopTime(cvode_mem, input_params_->tstop_);
  CVodeSetMaxOrd(cvode_mem, kBDFMaxOrder);
}

void VariableTimeStep::Clear(Branch *branch) {
  assert(branch->soma_);
  VariableTimeStep *vardt = (VariableTimeStep *)branch->interpolator_;
  vardt->~VariableTimeStep();
  branch->interpolator_ = nullptr;
}

void VariableTimeStep::PrintStatistics(const Branch *branch) {
  VariableTimeStep *vardt = (VariableTimeStep *)branch->interpolator_;
  CVodeMem cvode_mem = vardt->cvode_mem_;

  long num_steps = -1, num_rhs = -1, num_roots = -1, num_others = 0;
  CVodeGetNumSteps(cvode_mem, &num_steps);
  CVodeGetNumGEvals(cvode_mem, &num_roots);
  switch (input_params_->interpolator_) {
    case InterpolatorIds::kCvodePreConditionedDiagSolver:
      CVodeGetNumRhsEvals(cvode_mem, &num_rhs);
      CVodeGetNumNonlinSolvIters(cvode_mem, &num_others);
      printf(
          "-- Neuron %d completed. steps: %d, rhs: %d, pre-cond. solves: %d, "
          "roots: %d\n",
          branch->soma_->gid_, num_steps, num_rhs, num_others, num_roots);
      break;
    case InterpolatorIds::kCvodeDenseMatrix:
      CVDlsGetNumJacEvals(cvode_mem, &num_others);
      CVDlsGetNumRhsEvals(cvode_mem, &num_rhs);
      printf(
          "-- Neuron %d completed. steps: %d, rhs: %d, jacobians: %d, roots: "
          "%d\n",
          branch->soma_->gid_, num_steps, num_rhs, num_others, num_roots);
      break;
    case InterpolatorIds::kCvodeDiagonalMatrix:
      CVDiagGetNumRhsEvals(cvode_mem, &num_rhs);
      printf("-- Neuron %d completed. steps: %d, rhs: %d, roots: %d\n",
             branch->soma_->gid_, num_steps, num_rhs, num_roots);
      break;
    case InterpolatorIds::kCvodeSparseMatrix:
      CVDlsGetNumJacEvals(cvode_mem, &num_others);
      CVDlsGetNumRhsEvals(cvode_mem, &num_rhs);
      // CVSlsGetNumJacEvals(cvode_mem, &num_jacob_evals);
      // CVSlsGetNumRhsEvals(cvode_mem, &num_rhs_evals);
      break;
  }
}

hpx_t VariableTimeStep::StepTo(Branch *branch, const double tstop) {
  VariableTimeStep *vardt = (VariableTimeStep *)branch->interpolator_;
  NrnThread *nt = branch->nt_;
  CVodeMem cvode_mem = vardt->cvode_mem_;
  hpx_t spikes_lco = HPX_NULL;
  int roots_found[1];  // AP-threshold
  int flag = CV_ERR_FAILURE;
  floble_t event_group_ms = input_params_->cvode_event_group_ + 1e-12;

  double cvode_tstop = -1;
  while (nt->_t < tstop) {
    // delivers all events whithin the next delivery-time-window
    branch->DeliverEvents(nt->_t + event_group_ms);

    // get tout as time of next undelivered event (if any)
    hpx_lco_sema_p(branch->events_queue_mutex_);
    if (!branch->events_queue_.empty())
      cvode_tstop = std::min(tstop, branch->events_queue_.top().first);
    hpx_lco_sema_v_sync(branch->events_queue_mutex_);

    // call CVODE method: steps until reaching tout, or hitting root;
    while (nt->_t < cvode_tstop) {
      // perform several steps until hitting cvode_stop, or spiking
      flag = CVode(cvode_mem, cvode_tstop, vardt->y_, &(nt->_t), CV_NORMAL);

      // CVODE succeeded and roots found
      if (flag == CV_ROOT_RETURN) {
        flag = CVodeGetRootInfo(cvode_mem, roots_found);
        assert(flag == CV_SUCCESS);
        assert(roots_found[0] != 0);  // only root: AP threshold reached
        if (roots_found[0] > 0)       // AP-threshold reached from below (>0)
        {
          // if root found, integrator time is now at time of root
          hpx_t new_spikes_lco = branch->soma_->SendSpikes(nt->_t);
          if (new_spikes_lco) {
            // make sure only an AP occurred in between
            assert(!spikes_lco);
            spikes_lco = new_spikes_lco;
          }
        }
      } else {
        // error test failed repeatedly or with |h| = hmin.
        if (flag == CV_ERR_FAILURE) {
        }
        // convergence test failed too many times, or min-step was reached
        else if (flag == CV_CONV_FAILURE) {
        }
        // must have reached end, or we are into an unhandled error
        else {
          // success: nt->_t may have reached cvode_tstop or not;
        }
      }
    }
  }
  /* system deadlocks if neurons cant sent notifications a bit ahead.
   * They all depend on a step message that noone can send because its
   * to be send in the future. So now we send notifications for upcoming
   * messages */
  synchronizer_->StepSync(branch, 0);

  assert(fabs(nt->_t - tstop) < 0.001);  // equal
  return spikes_lco;
}
