#include "neurox/interpolators/variable_time_step.h"

#include <set>

using namespace neurox;
using namespace neurox::interpolators;
using namespace neurox::tools;

void VariableTimeStep::ScatterY(Branch *branch, N_Vector y) {
  const VariableTimeStep *vardt = (VariableTimeStep *)branch->vardt_;
  const double *y_data = NV_DATA_S(y);
  double **var_map = vardt->state_var_map_;
  for (int i = 0; i < vardt->equations_count_; i++) *(var_map[i]) = y_data[i];
}

void VariableTimeStep::GatherY(Branch *branch, N_Vector y) {
  const VariableTimeStep *vardt = (VariableTimeStep *)branch->vardt_;
  double *y_data = NV_DATA_S(y);
  double **var_map = vardt->state_var_map_;
  for (int i = 0; i < vardt->equations_count_; i++) y_data[i] = *(var_map[i]);
}

void VariableTimeStep::ScatterYdot(Branch *branch, N_Vector ydot) {
  const VariableTimeStep *vardt = (VariableTimeStep *)branch->vardt_;
  const double *ydot_data = NV_DATA_S(ydot);
  double **dv_map = vardt->state_dv_map_;
  for (int i = 0; i < vardt->equations_count_; i++) *(dv_map[i]) = ydot_data[i];
}

void VariableTimeStep::GatherYdot(Branch *branch, N_Vector ydot) {
  const VariableTimeStep *vardt = (VariableTimeStep *)branch->vardt_;
  double *ydot_data = NV_DATA_S(ydot);
  double **dv_map = vardt->state_dv_map_;
  for (int i = 0; i < vardt->equations_count_; i++) ydot_data[i] = *(dv_map[i]);
}

int VariableTimeStep::RootFunction(realtype t, N_Vector y, realtype *gout,
                                   void *user_data) {
  Branch *branch = (Branch *)user_data;

  // get offset and voltage of Axon Initial Segment
  int ais_offset = branch->thvar_ptr_ - branch->nt_->_actual_v;
  double v_ais = NV_Ith_S(y, ais_offset);

  // How it works: when gout[x] is zero, a root is found
  assert(v_ais >= -100 && v_ais < 50);
  gout[0] = v_ais - branch->soma_->threshold_;  // AP threshold reached
  return CV_SUCCESS;
}

int VariableTimeStep::RHSFunction(realtype t, N_Vector y, N_Vector ydot,
                                  void *user_data) {
  Branch *branch = (Branch *)user_data;
  VariableTimeStep *vardt = (VariableTimeStep *)branch->vardt_;
  NrnThread *nt = branch->nt_;
  realtype *ydot_data = NV_DATA_S(ydot);
  realtype *y_data = NV_DATA_S(y);

  //////// occvode.cpp: Cvode::fun_thread_transfer_part1

  const double h = vardt->cvode_mem_->cv_h;
  // nt->_dt = h==0 ? VariableTimeStep::kMinStepSize : h;
  nt->_dt = h == 0 ? 1e-8 : h;
  nt->cj = 1 / nt->_dt;
  nt->_t = t;

  // Updates internal states of continuous point processes (vecplay)
  // e.g. stimulus. vecplay->pd points to a read-only var used by
  // point proc mechanisms' nrn_current function
  // cvtrset.cpp :: CVode::fun_thread_transfer_part1()

  branch->FixedPlayContinuous(nt->_t);

  // copies V and state-vars from CVODES to NrnThread
  VariableTimeStep::ScatterY(branch, y);

  // for (int i=0; i<NV_LENGTH_S(vardt->y_); i++)
  //    fprintf(stderr, "== t=%.5f\ty[%d]=%.12f (init1)\n", t, i, y_data[i]);

  // start of occvode.cpp :: nocap_v
  solver::HinesSolver::ResetRHSandDNoCapacitors(branch, vardt->no_cap_);

  // sum mech-instance contributions to D and RHS on no-caps
  branch->CallModFunction(Mechanism::ModFunctions::kCurrent,
                          vardt->no_cap_->no_caps_ml_);  // rhs
  branch->CallModFunction(Mechanism::ModFunctions::kJacob,
                          vardt->no_cap_->no_caps_ml_);  // lhs

  solver::HinesSolver::SetupMatrixVoltageNoCapacitors(branch, vardt->no_cap_);

  // for (int i=0; i<NV_LENGTH_S(vardt->y_); i++)
  //    fprintf(stderr, "== t=%.5f\ty[%d]=%.12f (nocap_v)\n", t, i, y_data[i]);

  //////// ocvode2.cpp: Cvode::fun_thread_transfer_part2

  // fprintf(stderr, "PHASE2====\n");

  // cvtrset.cpp :: CVode::rhs
  solver::HinesSolver::ResetArray(branch, nt->_actual_rhs);

  // cvtrset.cpp :: CVode::rhs() -> rhs_memb(z.cv_memb_list_)
  // sum mech-instance contributions to D and RHS on caps
  branch->CallModFunction(Mechanism::ModFunctions::kCurrent);

  // add parent and children axial currents (A*dv and B*dv) to RHS
  // cvtrset.cpp :: CVode::rhs()
  solver::HinesSolver::SetupMatrixRHS(branch);

  // update mechanisms state (eg opening vars and derivatives)
  // cvtrset.cpp :: CVode::fun_thread_transfer_part2() -> do_ode()
  branch->CallModFunction(Mechanism::ModFunctions::kODESpec);

  // TODO missing logn term difus?
  // divide RHS by Cm and compute capacity current
  // cvtrset.cpp :: CVode::fun_thread_transfer_part2() -> nrn_div_capacity()
  branch->CallModFunction(Mechanism::ModFunctions::kDivCapacity);

  // copies dV/dt (RHS) and state-vars-derivative to CVODES
  VariableTimeStep::GatherYdot(branch, ydot);

  // for (int i=0; i<NV_LENGTH_S(vardt->y_); i++)
  //    fprintf(stderr, "== t=%.5f\tydot[%d]=%.12f (end2)\n", t, i,
  //    ydot_data[i]);

  return CV_SUCCESS;
}

// The solve function must solve the linear system Mx=b, where M is an
// approximation of the Newton matrix, I - gamma J, where J = df/dy,
// and the rhs-vector b is an input. Called once per newton, thus
// several times per time-step. (occvode.cpp: Cvode::solvex_thread)
int VariableTimeStep::NeuronLinearSolverFunction(CVodeMem cv_mem, N_Vector b,
                                                 N_Vector weight, N_Vector ycur,
                                                 N_Vector fcur) {
  // b is the right-hand-side vector, solution to be returned in b
  // ycur contains vector approximations to y(t_n)
  // ycur contains vector approximations to f(t_n, ycur)

  Branch *branch;
  NrnThread *nt = branch->nt_;
  nt->_dt = cv_mem->cv_gamma;
  nt->cj = 1 / dt;

  // Cvode::lhs()
  solver::HinesSolver::ResetArray(branch, nt->_actual_d);
  branch->CallModFunction(Mechanism::ModFunctions::kJacob);
  branch->CallModFunction(Mechanism::ModFunctions::kJacobCapacitance);
  nt->_actual_d[0] -= nt->_actual_b[0];
  solver::HinesSolver::SetupMatrixDiagonal(branch);
  // end of Cvode::lhs()

  ScatterYdot(branch, b);
  branch->CallModFunction(Mechanism::ModFunctions::kMulCapacity);
  solver::HinesSolver::ResetRHSNoCapacitors(branch, branch->vardt_);
  solver::HinesSolver::BackwardTriangulation(branch);
  solver::HinesSolver::ForwardSubstituion(branch);
  branch->CallModFunction(Mechanism::ModFunctions::kODEMatsol);
  GatherYdot(branch, b);
  return CV_SUCCESS;
}

// jacobian routine: compute J(t,y) = df/dy
int VariableTimeStep::JacobianDense(long int N, realtype t, N_Vector y,
                                    N_Vector fy, DlsMat J, void *user_data,
                                    N_Vector, N_Vector, N_Vector) {
  realtype **jac = J->cols;
  Branch *branch = (Branch *)user_data;
  VariableTimeStep *vardt = (VariableTimeStep *)branch->vardt_;
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
  // solver::HinesSolver::ResetMatrixV(branch);

  // does nothing so far (called before as part of nrn_current)
  // cvtrset.cpp :: CVode:: lhs -> lhs_memb()
  branch->CallModFunction(Mechanism::ModFunctions::kJacob);

  // D = D + cfac * .001 * _nt->cj // decay of D?
  branch->CallModFunction(Mechanism::ModFunctions::kJacobCapacitance);

  /*
  //add axial currents to D (d[i] -= b[i] and d[p[i]] -= a[i])
  solver::HinesSolver::SetupMatrixDiagonal(branch);

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
  assert(dv_offset == N);
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
      no_cap_(nullptr) {}

VariableTimeStep::~VariableTimeStep() {
  N_VDestroy_Serial(y_);
  CVodeFree((void **)(&cvode_mem_));
  delete no_cap_;
}

VariableTimeStep::NoCapacitor::~NoCapacitor() {
  delete child_ids_;
  delete node_ids_;
  // because these point to the same memory as
  // nt->data, we only delete main array
  delete[] no_caps_ml_;
  // Branch::DeleteMembList(memb_list_);
}

VariableTimeStep::NoCapacitor::NoCapacitor(const Branch *branch) {
  NrnThread *nt = branch->nt_;
  Memb_list *capac_instances = &branch->mechs_instances_[mechanisms_map_[CAP]];

  this->node_count_ = nt->end - capac_instances->nodecount;
  this->node_ids_ = new int[this->node_count_];
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
      this->node_ids_[no_cap_count++] = i;

    // if parent node is not a capacitors node
    if (i > 0 &&
        capacitor_ids.find(nt->_v_parent_index[i]) == capacitor_ids.end())
      child_ids.push_back(i);
  }

  // create childs ids and count of no-cap parents
  this->child_count_ = child_ids.size();
  this->child_ids_ = new int[this->child_count_];
  memcpy(this->child_ids_, child_ids.data(), child_ids.size() * sizeof(int));
  assert(this->node_count_ == no_cap_count);

  // occvode.cpp::new_no_cap_memb()
  input::DataLoader::GetMembListsOrderedByCapacitorsOrNot(
      branch, capacitor_ids, &(this->no_caps_ml_), nullptr);
}

// Neuron :: occvode.cpp :: init_global()
hpx_action_t VariableTimeStep::Init = 0;
int VariableTimeStep::Init_handler() {
  NEUROX_MEM_PIN(neurox::Branch);
  assert(local->vardt_ == nullptr);
  local->vardt_ = new VariableTimeStep();
  VariableTimeStep *vardt = (VariableTimeStep *)local->vardt_;
  CVodeMem &cvode_mem = vardt->cvode_mem_;
  int cap_count = local->mechs_instances_[mechanisms_map_[CAP]].nodecount;
  NrnThread *&nt = local->nt_;

  int flag = CV_ERR_FAILURE;

  // some methods from Branch::Finitialize
  // local->Finitialize2();
  local->CallModFunction(Mechanism::ModFunctions::kThreadTableCheck);
  local->InitVecPlayContinous();
  local->DeliverEvents(t);
  for (int n = 0; n < local->nt_->end; n++)
    local->nt_->_actual_v[n] = input_params_->voltage_;
  local->CallModFunction(Mechanism::ModFunctions::kBeforeInitialize);
  local->CallModFunction(Mechanism::ModFunctions::kInitialize);
  local->CallModFunction(Mechanism::ModFunctions::kAfterInitialize);
  local->CallModFunction(Mechanism::ModFunctions::kBeforeStep);
  local->DeliverEvents(t);

  // equations: capacitors + mechanisms * states
  int &equations_count = vardt->equations_count_;
  equations_count = cap_count;
  for (int m = 0; m < neurox::mechanisms_count_; m++) {
    Mechanism *mech = mechanisms_[m];
    // if (mech->type_==137 || mech->type_==139) continue; //TODO delete
    Memb_list *mech_instances = &local->mechs_instances_[m];
    equations_count += mech_instances->nodecount * mech->state_vars_->count_;
  }

  ///// create initial state y_0 for state array y

  // create map from y and dy to NrnThread->data
  vardt->state_var_map_ = new double *[equations_count]();
  vardt->state_dv_map_ = new double *[equations_count]();

  int var_offset = 0;
  Memb_list *capac_instances = &local->mechs_instances_[mechanisms_map_[CAP]];
  for (int c = 0; c < capac_instances->nodecount; c++) {
    int compartment_id = capac_instances->nodeindices[c];
    vardt->state_var_map_[var_offset] = &(nt->_actual_v[compartment_id]);
    vardt->state_dv_map_[var_offset] = &(nt->_actual_rhs[compartment_id]);
    var_offset++;
  }

  // collect information about non-capacitors nodes
  vardt->no_cap_ = new NoCapacitor(local);

  // build remaining map of state vars
  for (int m = 0; m < neurox::mechanisms_count_; m++) {
    Mechanism *mech = mechanisms_[m];
    Memb_list *mech_instances = &local->mechs_instances_[m];
    int ml_data_offset = 0;

// if (mech->type_==137 || mech->type_==139) continue; //TODO delete

/*
fprintf(stderr, "Mech %d , states %d*%d (neq=%d)\n",
       mech->type_, mech_instances->nodecount,
       mech->state_vars_->count_, var_offset);
*/

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
        int state_var_offset =
            ml_data_offset +
            Vectorizer::SizeOf(mech_instances->nodecount) * state_var_index + n;
        int state_dv_offset =
            ml_data_offset +
            Vectorizer::SizeOf(mech_instances->nodecount) * state_dv_index + n;
#endif
        // Note: SoA is OK but without padded memory on
        // CVODES maps, because we'd need to add empty
        // elements to CVODES y and y' arrays
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
  VariableTimeStep::GatherY(local, vardt->y_);

  // absolute tolerance array (low for voltages, high for mech states)
  vardt->absolute_tolerance_ = N_VNew_Serial(equations_count);
  for (int i = 0; i < equations_count; i++) {
    double tol = i < cap_count ? kAbsToleranceVoltage : kAbsToleranceMechStates;
    NV_Ith_S(vardt->absolute_tolerance_, i) = tol;
  }

  // CVodeCreate creates an internal memory block for a problem to
  // be solved by CVODES, with Backward Differentiation (or Adams)
  // and Newton solver (recommended for stiff problems, see header)
  cvode_mem = (CVodeMem)CVodeCreate(CV_BDF, CV_NEWTON);

  // from cvodeobj.cpp :: cvode_init()
  vardt->cvode_mem_->cv_gamma = 0.;
  vardt->cvode_mem_->cv_h = 0.;

  // CVodeInit allocates and initializes memory for a problem
  double t0 = input_params_->tstart_;
  flag = CVodeInit(cvode_mem, VariableTimeStep::RHSFunction, t0, vardt->y_);
  assert(flag == CV_SUCCESS);

  // specify integration tolerances. MUST be called before CVode.
  flag = CVodeSVtolerances(cvode_mem, kRelativeTolerance,
                           vardt->absolute_tolerance_);
  assert(flag == CV_SUCCESS);

  // specify user data to be used on functions as void* user_data_ptr;
  flag = CVodeSetUserData(cvode_mem, local);
  assert(flag == CV_SUCCESS);

  // specify root func. and roots (AP-threshold reached from below)
  int roots_direction[1] = {1};
  flag = CVodeRootInit(cvode_mem, 1, VariableTimeStep::RootFunction);
  CVodeSetRootDirection(cvode_mem, roots_direction);
  assert(flag == CV_SUCCESS);

  // initializes the memory record and sets various function
  // fields specific to the dense linear solver module.
  // Note: direct solvers give the solution (LU-decomposition, etc)
  // Indirect solvers require iterations (eg Jacobi method)

  switch (input_params_->interpolator_) {
    case Interpolators::kCvodePreCondNeuronSolver:
      // CVODES guide chapter 8: Providing Alternate Linear Solver
      // Modules: only lsolve function is mandatory
      // (non-used functions need to be set to null)
      cvode_mem->cv_linit = nullptr;
      cvode_mem->cv_lsetup = nullptr;
      cvode_mem->cv_lfree = nullptr;
      cvode_mem->cv_setupNonNull = FALSE;
      cvode_mem->cv_lsolve = NeuronLinearSolverFunction;
      break;
    case Interpolators::kCvodeDenseMatrix:
      flag = CVDense(cvode_mem, equations_count);
      break;
    case Interpolators::kCvodeDiagMatrix:
      flag = CVDiag(cvode_mem);
      break;
    case Interpolators::kCvodeSparseMatrix:
      // TODO
      assert(0);

      // Requires installation of Superlumt or KLU
      // flag = CVSlsSetSparseJacFn(cvode_mem, nullptr);
      // int nnz = equations_count * equations_count;
      // flag = CVKLU(cvode_mem, 1, equations_count, nnz);
      // flag = CVSuperLUMT(cvode_mem, 1, equations_count, nnz);

      // if not found, uses a dense matrix with sparse values
      flag = CVDense(cvode_mem, equations_count);
      flag = CVDlsSetDenseJacFn(cvode_mem, VariableTimeStep::JacobianDense);
      break;
  }
  assert(flag == CV_SUCCESS);

  CVodeSetInitStep(cvode_mem, kMinStepSize);
  CVodeSetMinStep(cvode_mem, kMinStepSize);
  CVodeSetMaxStep(cvode_mem,neurox::min_delay_steps_);
  CVodeSetStopTime(cvode_mem, input_params_->tstop_);
  CVodeSetMaxOrd(cvode_mem, kBDFMaxOrder);

  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t VariableTimeStep::Run = 0;
int VariableTimeStep::Run_handler() {
  NEUROX_MEM_PIN(neurox::Branch);
  assert(local->soma_);
  VariableTimeStep *vardt = (VariableTimeStep *)local->vardt_;
  CVodeMem cvode_mem = vardt->cvode_mem_;
  NrnThread *nt = local->nt_;

  int roots_found[1];  // AP-threshold
  int flag = CV_ERR_FAILURE;
  realtype tout = input_params_->tstop_;

  while (nt->_t < input_params_->tstop_) {
    // delivers all events whithin the next delivery-time-window
    local->DeliverEvents(nt->_t + VariableTimeStep::kEventsDeliveryTimeWindow);

    // get tout as time of next undelivered event (if any)
    hpx_lco_sema_p(local->events_queue_mutex_);
    if (!local->events_queue_.empty()) {
      tout = local->events_queue_.top().first;
    }
    hpx_lco_sema_v_sync(local->events_queue_mutex_);
    tout = std::min(input_params_->tstop_, tout);

    // call CVODE method: steps until reaching/passing tout;
    flag = CVode(cvode_mem, tout, vardt->y_, &(nt->_t), CV_NORMAL);

    printf("Neuron %d: t = %0.4f V=%f\n", nt->id, nt->_t, *(local->thvar_ptr_));

    if (flag == CV_ROOT_RETURN)  // CVODE succeeded and roots found
    {
      flag = CVodeGetRootInfo(cvode_mem, roots_found);
      assert(flag == CV_SUCCESS);
      assert(roots_found[0] != 0);  // only possible root: AP threshold reached
      // if root found, integrator time is now at time of root
      //(+1 value ascending, -1 valued descending)
      if (roots_found[0] > 0)  // AP-threshold reached from below
      {
        hpx_t spikes_lco_ = local->soma_->SendSpikes(nt->_t);
      }
    }
  }

#ifndef NDEBUG
  // Final statistics output:
  long num_steps = -1, num_rhs_evals = -1, num_jacob_evals = 0;
  CVodeGetNumSteps(cvode_mem, &num_steps);

  switch (input_params_->interpolator_) {
    case Interpolators::kCvodePreCondNeuronSolver:
      CVSpilsGetNumRhsEvals(cvode_mem, &num_rhs_evals);
      CVSpilsGetNumPrecSolves(cvode_mem, &num_jacob_evals);
      printf(
          "-- Neuron %d completed. steps: %d, rhs: %d, pre-cond. solves: %d\n",
          local->soma_->gid_, num_steps, num_rhs_evals, num_jacob_evals);
      break;
    case Interpolators::kCvodeDenseMatrix:
      CVDlsGetNumJacEvals(cvode_mem, &num_jacob_evals);
      CVDlsGetNumRhsEvals(cvode_mem, &num_rhs_evals);
      printf("-- Neuron %d completed. steps: %d, rhs: %d, jacobians: %d\n",
             local->soma_->gid_, num_steps, num_rhs_evals, num_jacob_evals);
      break;
    case Interpolators::kCvodeDiagMatrix:
      CVDiagGetNumRhsEvals(cvode_mem, &num_rhs_evals);
      printf("-- Neuron %d completed. steps: %d, rhs: %d\n", local->soma_->gid_,
             num_steps, num_rhs_evals);
      break;
    case Interpolators::kCvodeSparseMatrix:
      CVDlsGetNumJacEvals(cvode_mem, &num_jacob_evals);
      CVDlsGetNumRhsEvals(cvode_mem, &num_rhs_evals);
      // CVSlsGetNumJacEvals(cvode_mem, &num_jacob_evals);
      // CVSlsGetNumRhsEvals(cvode_mem, &num_rhs_evals);
      break;
  }
#endif
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t VariableTimeStep::Clear = 0;
int VariableTimeStep::Clear_handler() {
  NEUROX_MEM_PIN(neurox::Branch);
  assert(local->soma_);
  VariableTimeStep *vardt = (VariableTimeStep *)local->vardt_;
  vardt->~VariableTimeStep();
  return neurox::wrappers::MemoryUnpin(target);
}

void VariableTimeStep::RegisterHpxActions() {
  wrappers::RegisterZeroVarAction(VariableTimeStep::Init,
                                  VariableTimeStep::Init_handler);
  wrappers::RegisterZeroVarAction(VariableTimeStep::Run,
                                  VariableTimeStep::Run_handler);
  wrappers::RegisterZeroVarAction(VariableTimeStep::Clear,
                                  VariableTimeStep::Clear_handler);
}
