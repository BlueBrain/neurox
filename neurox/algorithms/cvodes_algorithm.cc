#include "neurox/algorithms/cvodes_algorithm.h"

#include "cvodes/cvodes_impl.h"
#include "cvodes/cvodes_diag.h"

using namespace neurox;
using namespace neurox::algorithms;
using namespace neurox::tools;

CvodesAlgorithm::CvodesAlgorithm() {}

CvodesAlgorithm::~CvodesAlgorithm() {}

const AlgorithmId CvodesAlgorithm::GetId() { return AlgorithmId::kCvodes; }

const char *CvodesAlgorithm::GetString() { return "CVODES"; }

void CvodesAlgorithm::ScatterY(Branch *branch, N_Vector y) {
  BranchCvodes *branch_cvodes =
      (BranchCvodes *)branch->soma_->algorithm_metadata_;
  double *y_data = NV_DATA_S(y);

  const int compartments_count = branch->nt_->end;
  memcpy(branch->nt_->_actual_v, y_data, sizeof(double) * compartments_count);

  double ** var_map = branch_cvodes->state_var_map_;
  for (int i = 0; i < branch_cvodes->equations_count_ - compartments_count; i++)
    *(var_map[i]) = y_data[compartments_count+i];
}

void CvodesAlgorithm::GatherYdot(Branch *branch, N_Vector ydot) {
    BranchCvodes *branch_cvodes =
        (BranchCvodes *)branch->soma_->algorithm_metadata_;
  double *ydot_data = NV_DATA_S(ydot);

  const int compartments_count = branch->nt_->end;
  memcpy(ydot_data, branch->nt_->_actual_rhs,
         sizeof(double) * compartments_count);

  double ** dv_map = branch_cvodes->state_dv_map_;
  for (int i = 0; i < branch_cvodes->equations_count_ - compartments_count; i++)
      ydot_data[compartments_count+i] = *(dv_map[i]);
}

void CvodesAlgorithm::ScatterYdot(Branch *branch, N_Vector ydot) {
    BranchCvodes *branch_cvodes =
        (BranchCvodes *)branch->soma_->algorithm_metadata_;
  double *ydot_data = NV_DATA_S(ydot);

  const int compartments_count = branch->nt_->end;
  memcpy(branch->nt_->_actual_rhs, ydot_data,
         sizeof(double) * compartments_count);

  double ** dv_map = branch_cvodes->state_dv_map_;
  for (int i = 0; i < branch_cvodes->equations_count_ - compartments_count; i++)
    *(dv_map[i]) = ydot_data[compartments_count+i];
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
      (BranchCvodes *)branch->soma_->algorithm_metadata_;
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
      (BranchCvodes *)branch->soma_->algorithm_metadata_;
  NrnThread *nt = branch->nt_;
  assert(t == nt->_t);

  //RHS provided the righ-hand side
  //We will now compute the LHS
  //Newton Matrix: M = I - gamma*J

  //in neuron we solve Px=b with P approximates Id-gamma*J
  //const double gamma = ((CVodeMem)branch_cvodes->cvodes_mem_)->cv_gamma;
  //nt->_dt = gamma;
  //nt->cj = 1/nt->_dt;

  //here we use only J
  nt->_dt = 1;
  nt->cj = 1;

  // cvtrset.cpp :: CVode:: lhs
  solver::HinesSolver::ResetMatrixV(branch);

  //does nothing so far (called before as part of kCurrent)
  // cvtrset.cpp :: CVode:: lhs -> lhs_memb()
  branch->CallModFunction(Mechanism::ModFunctions::kJacob);

  //call nrn_jacob_cap, which sums contributions to D
  branch->CallModFunction(Mechanism::ModFunctions::kJacobCapacitance);

  //add axial currents to D (d[i] -= b[i] and d[p[i]] -= a[i])
  solver::HinesSolver::SetupMatrixDiagonal(branch);

  //in neuron: b or ydot is the RHS input vector
  //TODO CvodesAlgorithm::ScatterYdot(branch, ydot);

  branch->CallModFunction(Mechanism::ModFunctions::kMulCapacity);

  solver::HinesSolver::BackwardTriangulation(branch);
  solver::HinesSolver::ForwardSubstituion(branch);

  //if mechs stiff
  //branch->CallModFunction(Mechanism::ModFunctions::kODEMatsol);

  //TODO CvodesAlgorithm::GatherYdot(branch, ydot);

  /// end of neuron occvode.cpp::solvex_thread()

  const int a_offset = nt->_actual_a - nt->_data;
  const int b_offset = nt->_actual_b - nt->_data;
  const int d_offset = nt->_actual_d - nt->_data;
  const int compartments_count = nt->end;

  const double *a = &nt->_data[a_offset];
  const double *b = &nt->_data[b_offset];
  const double *d = &nt->_data[d_offset];
  const int *p = nt->_v_parent_index;

  // Jacobian for main current equation:
  // d(dV/dt) / dV     = sum_i g_i x_i  + D  //mechs currents + D
  // d(dV/dt) / dV_p   = -A          //parent compartment
  // d(dV/dt) / dV_c_i = -B_c_i  //children_i compartment

  for (int n = 0; n < compartments_count; n++) {
    jac[n][n] = d[n];     // D = d (dV_n/dt) /dV_n
    if (n == 0) continue;
//    jac[p[n]][n] = a[n];  // A = d (dV_p/dt) /dV_n
//    jac[n][p[n]] = b[n];  // B = d (dV_n/dt) /dV_p
  }

  // get new derivative of mechs states (ions do not have,
  // capacitance added already to D in nrn_jacob_capacitance)
  // wrong: matsol does m = b/x for diagonal matrix solver
  // branch->CallModFunction(Mechanism::ModFunctions::kODEMatsol);

  //int compartment_id=-1;
  int dv_offset = compartments_count;
  for (int m = 0; m < neurox::mechanisms_count_; m++) {
    Mechanism *mech = mechanisms_[m];
    Memb_list *mech_instances = &branch->mechs_instances_[m];
    for (int n = 0; n < mech_instances->nodecount; n++)
      for (int s = 0; s < mech->state_vars_->count_; s++) {
        //Reminder: state vars Dm represents d(dm_i/dt) / dm_j
        jac[dv_offset][dv_offset] = *(branch_cvodes->state_dv_map_[dv_offset - compartments_count]);
        dv_offset++;
      }
  }
  assert(dv_offset==N);
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
      equations_count_(-1),
      state_var_map_(nullptr),
      state_dv_map_(nullptr),
      y_(nullptr),
      spikes_lco_(HPX_NULL)
{}

CvodesAlgorithm::BranchCvodes::~BranchCvodes() {
  N_VDestroy_Serial(y_);   /* Free y vector */
  CVodeFree(&cvodes_mem_); /* Free integrator memory */
  delete[] data_bak_;
}

//Neuron :: occvode.cpp :: init_global()
hpx_action_t CvodesAlgorithm::BranchCvodes::Init = 0;
int CvodesAlgorithm::BranchCvodes::Init_handler() {
  NEUROX_MEM_PIN(neurox::Branch);
  assert(local->soma_->algorithm_metadata_);
  BranchCvodes *branch_cvodes =
      (BranchCvodes*) local->soma_->algorithm_metadata_;
  void *&cvodes_mem = branch_cvodes->cvodes_mem_;
  int &equations_count = branch_cvodes->equations_count_;
  int compartments_count = local->nt_->end;
  NrnThread *&nt = local->nt_;

  int flag = CV_ERR_FAILURE;

  //used in all ode_matsol, nrn_state and nrn_init. Set at runtime by RHSFunction
  local->nt_->_dt = 1.0; //this shouldnt be necessary

  // calling same methods as Algorithm::FixedStepInit()
  local->Finitialize2();
  local->CallModFunction(Mechanism::ModFunctions::kThreadTableCheck);

  // equations: voltages per compartments + mechanisms * states
  equations_count = compartments_count;
  for (int m = 0; m < neurox::mechanisms_count_; m++)
  {
      Mechanism * mech=mechanisms_[m];
      Memb_list* mech_instances = &local->mechs_instances_[m];
      equations_count += mech_instances->nodecount * mech->state_vars_->count_;
  }

  ///// create initial state y_0 for state array y

  //set voltages first
  floble_t *y_data = new floble_t[equations_count];
  for (int i = 0; i < compartments_count; i++)
      y_data[i] = nt->_actual_v[i];

  // create map from y and dy to NrnThread->data (mech-states)
  // and set initial values for opening variables
  branch_cvodes->state_var_map_ =
      new double *[equations_count - compartments_count];
  branch_cvodes->state_dv_map_ =
      new double *[equations_count - compartments_count];

  int state_var_count = 0;
  for (int m = 0; m < neurox::mechanisms_count_; m++) {
    Mechanism *mech = mechanisms_[m];
    Memb_list *mech_instances = &local->mechs_instances_[m];
    int ml_data_offset = 0;
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
        //TODO this should be stored in SoA as well for LAYOUT=1
        assert(state_var_offset < Vectorizer::SizeOf(mech_instances->nodecount) * mech->data_size_);
        assert(state_dv_offset  < Vectorizer::SizeOf(mech_instances->nodecount) * mech->data_size_);
        branch_cvodes->state_var_map_[state_var_count] =
            &(mech_instances->data[state_var_offset]);
        branch_cvodes->state_dv_map_[state_var_count] =
            &(mech_instances->data[state_dv_offset]);
        y_data[compartments_count + state_var_count] =
            mech_instances->data[state_var_offset];
        state_var_count++;
      }
    }
    ml_data_offset += Vectorizer::SizeOf(mech_instances->nodecount) * mech->data_size_;
  }
  branch_cvodes->y_ = N_VMake_Serial(equations_count, y_data);
  assert(state_var_count == equations_count-compartments_count);

  // absolute tolerance array (low for voltages, high for mech states)
  branch_cvodes->absolute_tolerance_ = N_VNew_Serial(equations_count);
  for (int i = 0; i < equations_count; i++)
  {
    double tol = i<compartments_count ? kAbsToleranceVoltage : kAbsToleranceMechStates;
    NV_Ith_S(branch_cvodes->absolute_tolerance_, i) = tol;
  }

  // CVodeCreate creates an internal memory block for a problem to
  // be solved by CVODES, with Backward Differentiation (or Adams)
  // and Newton solver (recommended for stiff problems, see header)
  cvodes_mem = CVodeCreate(CV_BDF, CV_NEWTON);

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
  branch_cvodes->data_bak_ = new double[local->nt_->_ndata];
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
  flag = CVDlsSetDenseJacFn(cvodes_mem, CvodesAlgorithm::JacobianDense);
#else  // sparse colver
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

  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t CvodesAlgorithm::BranchCvodes::Run = 0;
int CvodesAlgorithm::BranchCvodes::Run_handler() {
  NEUROX_MEM_PIN(neurox::Branch);
  assert(local->soma_);
  BranchCvodes *branch_cvodes =
      (BranchCvodes *)local->soma_->algorithm_metadata_;
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
    /*
    N_Vector estimated_local_errors =
    N_VNewEmpty_Serial(branch_cvodes->equations_count_);
    CVodeGetEstLocalErrors(cvodes_mem, estimated_local_errors);
    */
  }

  // Final statistics output:
  long num_steps = -1, num_jacob_evals = -1, num_rhs_evals = -1;
  realtype last_step_size = -1;

  CVodeGetNumSteps(cvodes_mem, &num_steps);
  CVodeGetLastStep(cvodes_mem, &last_step_size);
#if NEUROX_CVODES_JACOBIAN_SOLVER == 0
  CVDlsGetNumJacEvals(cvodes_mem, &num_jacob_evals);
  CVDlsGetNumRhsEvals(cvodes_mem, &num_rhs_evals);
#else
  CVSlsGetNumJacEvals(cvodes_mem, &num_jacob_evals);
  CVSlsGetNumRhsEvals(cvodes_mem, &num_rhs_evals);
#endif

  printf("- num_steps: %d\n", num_steps);
  printf("- num_jacob_evals: %d\n", num_jacob_evals);
  printf("- num_rhs_evals: %d\n", num_rhs_evals);
  printf("- last_step_size: %d\n", last_step_size);

  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t CvodesAlgorithm::BranchCvodes::Clear = 0;
int CvodesAlgorithm::BranchCvodes::Clear_handler() {
  NEUROX_MEM_PIN(neurox::Branch);
  assert(local->soma_);
  BranchCvodes *branch_cvodes =
      (BranchCvodes *)local->soma_->algorithm_metadata_;
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
