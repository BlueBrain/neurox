#include "neurox/algorithms/cvodes_algorithm.h"

#include "cvodes/cvodes_impl.h"
#include "cvodes/cvodes_diag.h"

#include <set>

using namespace neurox;
using namespace neurox::algorithms;
using namespace neurox::tools;

CvodesAlgorithm::CvodesAlgorithm() {}

CvodesAlgorithm::~CvodesAlgorithm() {}

const SyncAlgorithms CvodesAlgorithm::GetId() { return SyncAlgorithms::kAllReduce; }

const char *CvodesAlgorithm::GetString() { return "CVODES"; }

void CvodesAlgorithm::ScatterY(Branch *branch, N_Vector y) {
  BranchCvodes *branch_cvodes =
      (BranchCvodes *)branch->branch_cvodes_;
  const double *y_data = NV_DATA_S(y);
  double ** var_map = branch_cvodes->state_var_map_;
  for (int i = 0; i < branch_cvodes->equations_count_; i++)
    *(var_map[i]) = y_data[i];
}

void CvodesAlgorithm::GatherY(Branch *branch, N_Vector y) {
  BranchCvodes *branch_cvodes =
        (BranchCvodes *)branch->branch_cvodes_;
  double *y_data = NV_DATA_S(y);
  double ** var_map = branch_cvodes->state_var_map_;
  for (int i = 0; i < branch_cvodes->equations_count_; i++)
      y_data[i] = *(var_map[i]);
}

void CvodesAlgorithm::ScatterYdot(Branch *branch, N_Vector ydot) {
  BranchCvodes *branch_cvodes =
        (BranchCvodes *)branch->branch_cvodes_;
  const double *ydot_data = NV_DATA_S(ydot);
  double ** dv_map = branch_cvodes->state_dv_map_;
  for (int i = 0; i < branch_cvodes->equations_count_; i++)
    *(dv_map[i]) = ydot_data[i];
}

void CvodesAlgorithm::GatherYdot(Branch *branch, N_Vector ydot) {

  //on initialization we call RHS with ydot==NULL
  if (ydot==nullptr) return;

  BranchCvodes *branch_cvodes =
      (BranchCvodes *)branch->branch_cvodes_;
  double *ydot_data = NV_DATA_S(ydot);
  double ** dv_map = branch_cvodes->state_dv_map_;
  for (int i = 0; i < branch_cvodes->equations_count_; i++)
      ydot_data[i] = *(dv_map[i]);
}

/// g root function to compute g_i(t,y) .
int CvodesAlgorithm::RootFunction(realtype t, N_Vector y, realtype *gout,
                                  void *user_data) {
  Branch *branch = (Branch*)user_data;

  //get offset and voltage of Axon Initial Segment
  int ais_offset = branch->thvar_ptr_ - branch->nt_->_actual_v;
  double v_ais = NV_Ith_S(y, ais_offset);

  // How it works: when gout[x] is zero, a root is found
  assert(v_ais >= -100 && v_ais < 30);
  gout[0] = v_ais - branch->soma_->threshold_;  // AP threshold reached
  return CV_SUCCESS;
}

/// f routine to compute ydot=f(t,y), i.e. compute new values of nt->data
/// from neuron::occvode.cpp::fun_thread(...)
int CvodesAlgorithm::RHSFunction(realtype t, N_Vector y, N_Vector ydot,
                                 void *user_data) {
  Branch *branch = (Branch*) user_data;
  BranchCvodes *branch_cvodes =
      (BranchCvodes *)branch->branch_cvodes_;
  NrnThread *nt = branch->nt_;
  realtype *ydot_data = NV_DATA_S(ydot);
  realtype *y_data = NV_DATA_S(y);

  //////// occvode.cpp: fun_thread_transfer_part1 /////////

  const double h = ((CVodeMem)branch_cvodes->cvodes_mem_)->cv_h;
  nt->_dt = h==0 ? CvodesAlgorithm::kMinStepSize : h;
  nt->cj = 1/nt->_dt;
  nt->_t = t;

  // Updates internal states of continuous point processes (vecplay)
  // e.g. stimulus. vecplay->pd points to a read-only var used by
  // point proc mechanisms' nrn_current function
  //cvtrset.cpp :: CVode::fun_thread_transfer_part1()
  branch->FixedPlayContinuous(nt->_t);

  //copies V and state-vars from CVODES to NrnThread
  CvodesAlgorithm::ScatterY(branch, y);

  double * yy_data = NV_DATA_S(branch_cvodes->y_);
  for (int i=0; i<NV_LENGTH_S(branch_cvodes->y_); i++)
      if (yy_data[i]!=0)
        fprintf(stderr, "y0[%d]=%.12f\n", i, yy_data[i]);

  //start of occvode.cpp :: nocap_v
  for (int i=0; i<branch_cvodes->no_cap_count; ++i)
  {
      int no_cap_id = branch_cvodes->no_cap_node[i];
      nt->_actual_d[i] = 0;
      nt->_actual_rhs[i] = 0;
  }

  //TODO call only on non-capacitance nodes!
  solver::HinesSolver::ResetNoCapacitanceRHSandD(branch, branch_cvodes);
  branch->CallModFunction(Mechanism::ModFunctions::kCurrent);
  branch->CallModFunction(Mechanism::ModFunctions::kJacob);
  solver::HinesSolver::NoCapacitanceVoltage(branch, branch_cvodes);

  //////// ocvode2.cpp: fun_thread_transfer_part2 ///////

  //cvtrset.cpp :: CVode::rhs
  solver::HinesSolver::ResetMatrixRHSandD(branch);

  //cvtrset.cpp :: CVode::rhs() -> rhs_memb()
  //sum mech-instance contributions to D and RHS
  branch->CallModFunction(Mechanism::ModFunctions::kCurrent);

  // add parent and children axial currents (A*dv and B*dv) to RHS
  // cvtrset.cpp :: CVode::rhs()
  solver::HinesSolver::SetupMatrixRHS(branch);

  //update mechanisms state (eg opening vars and derivatives)
  // cvtrset.cpp :: CVode::fun_thread_transfer_part2() -> do_ode()
  branch->CallModFunction(Mechanism::ModFunctions::kODESpec);

  // divide RHS by Cm and compute capacity current
  // cvtrset.cpp :: CVode::fun_thread_transfer_part2() -> nrn_div_capacity()
  branch->CallModFunction(Mechanism::ModFunctions::kDivCapacity);

  //copies dV/dt (RHS) and state-vars-derivative to CVODES
  CvodesAlgorithm::GatherYdot(branch, ydot);

  printf ("RHS neuron %d, t%.10f, V=%.10f\n",
          branch->soma_->gid_, t, *branch->thvar_ptr_);

  return CV_SUCCESS;
}

// jacobian routine: compute J(t,y) = df/dy
int CvodesAlgorithm::JacobianDense(long int N, realtype t, N_Vector y,
                                      N_Vector fy, DlsMat J,
                                      void *user_data, N_Vector, N_Vector,
                                      N_Vector) {
  realtype **jac = J->cols;
  Branch *branch = (Branch*) user_data;
  BranchCvodes *branch_cvodes =
      (BranchCvodes *)branch->branch_cvodes_;
  NrnThread *nt = branch->nt_;
  assert(t == nt->_t);

  //RHS provided the righ-hand side
  //We will now compute the LHS
  //Newton Matrix: M = I - gamma*J

  //in neuron we solve Px=b with P approximates Id-gamma*J
  //const double gamma = ((CVodeMem)branch_cvodes->cvodes_mem_)->cv_gamma;
  //nt->_dt = gamma;
  //nt->cj = 1/nt->_dt;

  //used in kJacobCapacitance
  nt->_dt = 1;
  nt->cj = 1;

  // if I reset V, then nrn_current is wrong!
  // cvtrset.cpp :: CVode:: lhs
  // solver::HinesSolver::ResetMatrixV(branch);

  //does nothing so far (called before as part of nrn_current)
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

  const int capacitances_count = nt->end;
  const double *a = nt->_actual_a;
  const double *b = nt->_actual_b;
  const double *d = nt->_actual_d;
  const int *p = nt->_v_parent_index;

  // Jacobian for main current equation:
  // d(dV/dt) / dV     = sum_i g_i x_i  + D  //mechs currents + D
  // d(dV/dt) / dV_p   = -A          //parent compartment
  // d(dV/dt) / dV_c_i = -B_c_i  //children_i compartment

  for (int n = 0; n < capacitances_count; n++) {
    jac[n][n] = d[n];     // D = d (dV_n/dt) /dV_n
    if (n == 0) continue;
    jac[p[n]][n] = a[n];  // A = d (dV_p/dt) /dV_n
    jac[n][p[n]] = b[n];  // B = d (dV_n/dt) /dV_p
  }

  // Matsol  Matsol solves (1 + dt*jacobian)*x = b
  // for direct diagonal matrix solver
  // branch->CallModFunction(Mechanism::ModFunctions::kODEMatsol);

  //int compartment_id=-1;
  int dv_offset = capacitances_count;
  for (int m = 0; m < neurox::mechanisms_count_; m++) {
    Mechanism *mech = mechanisms_[m];
    Memb_list *mech_instances = &branch->mechs_instances_[m];
    for (int n = 0; n < mech_instances->nodecount; n++)
      for (int s = 0; s < mech->state_vars_->count_; s++) {
        //Reminder: state vars Dm represents d(dm_i/dt) / dm_j
        jac[dv_offset][dv_offset] = *(branch_cvodes->state_dv_map_[dv_offset - capacitances_count]);
        dv_offset++;
      }
  }
  assert(dv_offset==N);
  printf ("Jac neuron %d, t%.10f\n", branch->soma_->gid_, t);
  double* y_data = NV_DATA_S(y);
  //for (int i=0; i<N; i++)
  //    printf("t=%f\t RHS[%d] = %.10f\t Jac[%d] = %.10f\n", t, i, y_data[i], i, jac[i][i]);

  return CV_SUCCESS;
}

/////////////////////// Algorithm abstract class ////////////////////////

void CvodesAlgorithm::Init() {
  neurox::wrappers::CallAllNeurons(CvodesAlgorithm::BranchCvodes::Init);
}

void CvodesAlgorithm::Clear() {
  neurox::wrappers::CallAllNeurons(CvodesAlgorithm::BranchCvodes::Clear);
}

double CvodesAlgorithm::Launch() {
  hpx_time_t now = hpx_time_now();
  neurox::wrappers::CallAllNeurons(CvodesAlgorithm::BranchCvodes::Run);
  return hpx_time_elapsed_ms(now) / 1e3;
}

void CvodesAlgorithm::StepBegin(Branch *) {}

void CvodesAlgorithm::StepEnd(Branch *b, hpx_t spikesLco) {}

void CvodesAlgorithm::Run(Branch *b, const void *args) {}

hpx_t CvodesAlgorithm::SendSpikes(Neuron *n, double tt, double) {
  return Neuron::SendSpikesAsync(n, tt);
}

//////////////////////////// BranchCvodes /////////////////////////

CvodesAlgorithm::BranchCvodes::BranchCvodes()
    : cvodes_mem_(nullptr),
      capacitances_count_(-1),
      equations_count_(-1),
      state_var_map_(nullptr),
      state_dv_map_(nullptr),
      y_(nullptr),
      spikes_lco_(HPX_NULL)
{}

CvodesAlgorithm::BranchCvodes::~BranchCvodes() {
  N_VDestroy_Serial(y_);   /* Free y vector */
  CVodeFree(&cvodes_mem_); /* Free integrator memory */
}

//Neuron :: occvode.cpp :: init_global()
hpx_action_t CvodesAlgorithm::BranchCvodes::Init = 0;
int CvodesAlgorithm::BranchCvodes::Init_handler() {
  NEUROX_MEM_PIN(neurox::Branch);
  assert(local->branch_cvodes_==nullptr);
  local->branch_cvodes_ = new algorithms::CvodesAlgorithm::BranchCvodes();
  BranchCvodes *branch_cvodes =
      (BranchCvodes*) local->branch_cvodes_;
  void *&cvodes_mem = branch_cvodes->cvodes_mem_;
  branch_cvodes->capacitances_count_ = local->mechs_instances_[mechanisms_map_[CAP]].nodecount;
  NrnThread *&nt = local->nt_;

  int flag = CV_ERR_FAILURE;

  // some methods from Branch::Finitialize
  local->CallModFunction(Mechanism::ModFunctions::kThreadTableCheck);
  local->InitVecPlayContinous();
  local->CallModFunction(Mechanism::ModFunctions::kThreadTableCheck);

  // equations: capacitances + mechanisms * states
  int & equations_count = branch_cvodes->equations_count_;
  equations_count =  branch_cvodes->capacitances_count_;
  for (int m = 0; m < neurox::mechanisms_count_; m++)
  {
      Mechanism * mech=mechanisms_[m];
      if (mech->type_==137 || mech->type_==139) continue; //TODO delete
      Memb_list* mech_instances = &local->mechs_instances_[m];
      equations_count += mech_instances->nodecount * mech->state_vars_->count_;
  }

  ///// create initial state y_0 for state array y

  // create map from y and dy to NrnThread->data
  branch_cvodes->state_var_map_ = new double *[equations_count]();
  branch_cvodes->state_dv_map_  = new double *[equations_count]();

  int var_offset = 0;
  Memb_list *capac_instances = &local->mechs_instances_[mechanisms_map_[CAP]];
  std::set<int> capacitance_ids;
  for (int c=0; c<capac_instances->nodecount; c++)
  {
      int compartment_id = capac_instances->nodeindices[c];
      capacitance_ids.insert(compartment_id);
      branch_cvodes->state_var_map_[var_offset] =
          &(nt->_actual_v[compartment_id]);
      branch_cvodes->state_dv_map_[var_offset] =
          &(nt->_actual_rhs[compartment_id]);
      var_offset++;
  }

  /* MEM corruption
  //based on previous RHS marks, build tree of nodes and children
  branch_cvodes->no_cap_count = nt->end - branch_cvodes->capacitances_count_;
  branch_cvodes->no_cap_child = new int[branch_cvodes->no_cap_count];
  branch_cvodes->no_cap_node  = new int[branch_cvodes->no_cap_count];
  int no_cap_offset=0;
  int no_cap_child_count=0;
  for (int i=0; i<nt->end; i++)
  {
      //if not a capacitance node
      if (capacitance_ids.find(i)!=capacitance_ids.end())
          branch_cvodes->no_cap_node[no_cap_offset++] = i;
      else //build list of non-caps
          branch_cvodes->no_cap_list.push_back(i);

      //if parent is not a capacitance node
      if (i > 0 && capacitance_ids.find(nt->_v_parent_index[i])!=capacitance_ids.end())
          branch_cvodes->no_cap_child[no_cap_child_count++] = i;
  }
  branch_cvodes->no_cap_child_count = no_cap_child_count;
  */

  //build remaining map of state vars
  for (int m = 0; m < neurox::mechanisms_count_; m++) {
    Mechanism *mech = mechanisms_[m];
    Memb_list *mech_instances = &local->mechs_instances_[m];
    int ml_data_offset = 0;

      if (mech->type_==137 || mech->type_==139) continue; //TODO delete

      fprintf(stderr, "Mech %d , states %d*%d (neq=%d)\n",
             mech->type_,
             mech_instances->nodecount, mech->state_vars_->count_,
             var_offset);
    for (int n = 0; n < mech_instances->nodecount; n++) {
      for (int s = 0; s < mech->state_vars_->count_; s++) {
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
            Vectorizer::SizeOf(mech_instances->nodecount) * state_var_index +n;
        int state_dv_offset =
            ml_data_offset +
            Vectorizer::SizeOf(mech_instances->nodecount) * state_dv_index +n;
#endif
        //TODO this should be stored in SoA as well for LAYOUT=1 ??
        assert(state_var_offset < Vectorizer::SizeOf(mech_instances->nodecount) * mech->data_size_);
        assert(state_dv_offset  < Vectorizer::SizeOf(mech_instances->nodecount) * mech->data_size_);
        branch_cvodes->state_var_map_[var_offset] =
            &(mech_instances->data[state_var_offset]);
        branch_cvodes->state_dv_map_[var_offset] =
            &(mech_instances->data[state_dv_offset]);
        var_offset++;
      }
    }
    ml_data_offset += Vectorizer::SizeOf(mech_instances->nodecount) * mech->data_size_;
  }
  assert(var_offset == equations_count);
  branch_cvodes->y_ = N_VNew_Serial(equations_count);
  CvodesAlgorithm::GatherY(local, branch_cvodes->y_);

  // absolute tolerance array (low for voltages, high for mech states)
  branch_cvodes->absolute_tolerance_ = N_VNew_Serial(equations_count);
  for (int i = 0; i < equations_count; i++)
  {
    double tol = i<branch_cvodes->capacitances_count_ ? kAbsToleranceVoltage : kAbsToleranceMechStates;
    NV_Ith_S(branch_cvodes->absolute_tolerance_, i) = tol;
  }

  // CVodeCreate creates an internal memory block for a problem to
  // be solved by CVODES, with Backward Differentiation (or Adams)
  // and Newton solver (recommended for stiff problems, see header)
  cvodes_mem = CVodeCreate(CV_BDF, CV_NEWTON);

  // from cvodeobj.cpp :: cvode_init()
  ((CVodeMem)branch_cvodes->cvodes_mem_)->cv_gamma = 0.;
  ((CVodeMem)branch_cvodes->cvodes_mem_)->cv_h = 0.;

  // CVodeInit allocates and initializes memory for a problem
  // In Neuron RHSFn is cvodeobj.cpp :: f_lvardt
  double t0 = input_params_->tstart_;
  flag = CVodeInit(cvodes_mem, CvodesAlgorithm::RHSFunction, t0, branch_cvodes->y_);
  assert(flag == CV_SUCCESS);

  // specify integration tolerances. MUST be called before CVode.
  flag = CVodeSVtolerances(cvodes_mem, kRelativeTolerance,
                           branch_cvodes->absolute_tolerance_);
  assert(flag == CV_SUCCESS);

  // specify user data to be used on functions as void* user_data_ptr;
  flag = CVodeSetUserData(cvodes_mem, local);
  assert(flag == CV_SUCCESS);

  // specify g as root function and roots
  int roots_direction[1] = {1};  // AP threshold reached for increasing voltage
  flag = CVodeRootInit(cvodes_mem, 1, CvodesAlgorithm::RootFunction);
  CVodeSetRootDirection(cvodes_mem, roots_direction);
  assert(flag == CV_SUCCESS);

// initializes the memory record and sets various function
// fields specific to the dense linear solver module.
// Note: direct solvers give the solution (LU-decomposition, etc)
// Indirect solvers require iterations (eg Jacobi method)

#if NEUROX_CVODES_JACOBIAN_SOLVER == 0  // approx diag Jac
  //no Jacobian provided, uses an approx. diagonal Jacobian
  flag = CVDiag(cvodes_mem);
#elif NEUROX_CVODES_JACOBIAN_SOLVER ==1   // dense colver
  flag = CVDense(cvodes_mem, equations_count);
  //TODO
  //flag = CVDlsSetDenseJacFn(cvodes_mem, CvodesAlgorithm::JacobianDense);
#else  // sparse solver
  // Requires installation of Superlumt or KLU
  flag =
      CVSlsSetSparseJacFn(cvode_mem_, CvodesAlgorithm::JacobianSparseFunction);
  int nnz = equations_count * equations_count;
#if NEUROX_CVODES_JACOBIAN_SOLVER == 2
  flag = CVKLU(cvode_mem, 1, equations_count, nnz);
#else  //==3
  flag = CVSuperLUMT(cvode_mem, 1, equations_count, nnz);
#endif
#endif

  assert(flag == CV_SUCCESS);

  //CVodeSetInitStep(cvodes_mem, kMinStepSize);
  CVodeSetMinStep(cvodes_mem, kMinStepSize);
  CVodeSetMaxStep(cvodes_mem,
                  CoreneuronAlgorithm::CommunicationBarrier::kCommStepSize);
  CVodeSetStopTime(cvodes_mem, input_params_->tstop_);
  CVodeSetMaxOrd(cvodes_mem, kBDFMaxOrder);

  // call one RHS (cvodeobj.cpp :: Cvode::cvode_init)
  GatherY(branch, branch_cvodes->y_);
  RHSFunction(input_params_->tstart_,
              branch_cvodes->y_,
              NULL, local);
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t CvodesAlgorithm::BranchCvodes::Run = 0;
int CvodesAlgorithm::BranchCvodes::Run_handler() {
  NEUROX_MEM_PIN(neurox::Branch);
  assert(local->soma_);
  BranchCvodes *branch_cvodes =
      (BranchCvodes *)local->branch_cvodes_;
  void *cvodes_mem = branch_cvodes->cvodes_mem_;
  NrnThread *nt = local->nt_;

  int roots_found[1];  // AP-threshold
  int flag = CV_ERR_FAILURE;
  realtype tout = input_params_->tstop_;

  while (nt->_t < input_params_->tstop_) {
    // delivers all events whithin the next delivery-time-window
    local->DeliverEvents(nt->_t + CvodesAlgorithm::kEventsDeliveryTimeWindow);

    // get tout as time of next undelivered event (if any)
    hpx_lco_sema_p(local->events_queue_mutex_);
    if (!local->events_queue_.empty()) {
        tout = local->events_queue_.top().first;
    }
    hpx_lco_sema_v_sync(local->events_queue_mutex_);
    tout = std::min(input_params_->tstop_, tout);

    // call CVODE method: steps until reaching/passing tout;
    flag = CVode(cvodes_mem, tout, branch_cvodes->y_, &(nt->_t), CV_NORMAL);

    double *v = NV_DATA_S(branch_cvodes->y_);
    printf("Neuron %d: t = %0.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", nt->id,
           nt->_t, v[0], v[1], v[2], v[3], v[4], v[5]);

    if (flag == CV_ROOT_RETURN)  // CVODE succeeded and roots found
    {
      flag = CVodeGetRootInfo(cvodes_mem, roots_found);
      assert(flag == CV_SUCCESS);
      assert(roots_found[0] != 0);  // only possible root: AP threshold reached
      // NOTE: if root is found, integrator time is now at time of root!
      //(+1 value ascending, -1 valued descending)
      if (roots_found[0] > 0)  // AP threshold reached from below
      {
        // TODO we can use a root to stop progress in time and wait for deliver?
        // use several hpx-all reduce to mark performance?
        branch_cvodes->spikes_lco_ = local->soma_->SendSpikes(nt->_t);
      }
    }
  }

  // Final statistics output:
  long num_steps = -1, num_rhs_evals = -1;
  CVodeGetNumSteps(cvodes_mem, &num_steps);
  //CVodeGetLastStep(cvodes_mem, &last_step_size);
#if NEUROX_CVODES_JACOBIAN_SOLVER == 0 //CVdiag (no jacobian)
  CVDiagGetNumRhsEvals(cvodes_mem, &num_rhs_evals);
  printf("- num_steps: %d, num_rhs_evals: %d\n",
         num_steps, num_rhs_evals);
#else
  long num_jacob_evals = -1;
  #if NEUROX_CVODES_JACOBIAN_SOLVER == 1
    CVDlsGetNumJacEvals(cvodes_mem, &num_jacob_evals);
    CVDlsGetNumRhsEvals(cvodes_mem, &num_rhs_evals);
  #else
    CVSlsGetNumJacEvals(cvodes_mem, &num_jacob_evals);
    CVSlsGetNumRhsEvals(cvodes_mem, &num_rhs_evals);
  #endif
  printf("- num_steps: %d,  num_jacob_evals: %d, num_rhs_evals: %d\n",
         num_steps, num_jacob_evals, num_rhs_evals);
#endif
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t CvodesAlgorithm::BranchCvodes::Clear = 0;
int CvodesAlgorithm::BranchCvodes::Clear_handler() {
  NEUROX_MEM_PIN(neurox::Branch);
  assert(local->soma_);
  BranchCvodes *branch_cvodes =
      (BranchCvodes *)local->branch_cvodes_;
  branch_cvodes->~BranchCvodes();
  return neurox::wrappers::MemoryUnpin(target);
}

void CvodesAlgorithm::BranchCvodes::RegisterHpxActions() {
  wrappers::RegisterZeroVarAction(CvodesAlgorithm::BranchCvodes::Init,
                                  CvodesAlgorithm::BranchCvodes::Init_handler);
  wrappers::RegisterZeroVarAction(CvodesAlgorithm::BranchCvodes::Run,
                                  CvodesAlgorithm::BranchCvodes::Run_handler);
  wrappers::RegisterZeroVarAction(CvodesAlgorithm::BranchCvodes::Clear,
                                  CvodesAlgorithm::BranchCvodes::Clear_handler);
}
