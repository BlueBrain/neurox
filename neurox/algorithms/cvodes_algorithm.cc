#include "neurox/algorithms/cvodes_algorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

CvodesAlgorithm::CvodesAlgorithm() {}

CvodesAlgorithm::~CvodesAlgorithm() {}

const AlgorithmId CvodesAlgorithm::GetId() { return AlgorithmId::kCvodes; }

const char *CvodesAlgorithm::GetString() { return "CVODES"; }

void CvodesAlgorithm::CopyYToVoltage(N_Vector y, Branch *branch) {
  BranchCvodes *branch_cvodes =
      (BranchCvodes *)branch->soma_->algorithm_metadata_;
  double *y_data = NV_DATA_S(y);

  int compartments_count = branch->nt_->end;
  memcpy(branch->nt_->_actual_v, y_data, sizeof(double) * compartments_count);

  double ** var_map = branch_cvodes->state_var_map_;
  for (int i = 0; i < branch_cvodes->equations_count_ - compartments_count; i++)
    *(var_map[i]) = y_data[i];
}

void CvodesAlgorithm::CopyRHSToYdot(Branch *branch, N_Vector ydot) {
    BranchCvodes *branch_cvodes =
        (BranchCvodes *)branch->soma_->algorithm_metadata_;
  double *ydot_data = NV_DATA_S(ydot);

  int compartments_count = branch->nt_->end;
  memcpy(ydot_data, branch->nt_->_actual_rhs,
         sizeof(double) * compartments_count);

  double ** dv_map = branch_cvodes->state_dv_map_;
  for (int i = 0; i < branch_cvodes->equations_count_ - compartments_count; i++)
      ydot_data[i] = *(dv_map[i]);
}

/// g root function to compute g_i(t,y) .
int CvodesAlgorithm::RootFunction(realtype t, N_Vector y, realtype *gout,
                                  void *user_data_ptr) {
  Branch *branch = ((BranchCvodes::UserData *)user_data_ptr)->branch_;

  //get offset and voltage of Axon Initial Segment
  int ais_offset = branch->thvar_ptr_ - branch->nt_->_actual_v;
  double v_ais = NV_Ith_S(y, ais_offset);

  // How it works: when gout[x] is zero, a root is found
  assert(v_ais >= -100 && v_ais < 30);
  gout[0] = v_ais - branch->soma_->threshold_;  // AP threshold reached
  return CV_SUCCESS;
}

/// f routine to compute ydot=f(t,y), i.e. compute new values of nt->data
int CvodesAlgorithm::RHSFunction(realtype t, N_Vector y, N_Vector ydot,
                                 void *user_data_ptr) {
  BranchCvodes::UserData *user_data = (BranchCvodes::UserData *)user_data_ptr;
  Branch *branch = user_data->branch_;
  NrnThread *nt = branch->nt_;
  realtype *ydot_data = NV_DATA_S(ydot);
  realtype *y_data = NV_DATA_S(y);

  nt->_t = t;   // nt->t used in FixedPlayContinuous
  nt->_dt = t - user_data->rhs_last_time_; // nt->dt used in ode_matsol1
  nt->cj = (input_params_->second_order_ ? 2.0 : 1.0) / nt->_dt;
  user_data->rhs_second_last_time_ = user_data->rhs_last_time_;
  user_data->rhs_last_time_ = t;

  // update vars in NrnThread->data described by our CVODES state
  CvodesAlgorithm::CopyYToVoltage(y, branch);

  /////////   Get new RHS and D from current state ///////////
  // Note: coreneuron computes RHS and jacobian (D) simultaneously.

  // Updates internal states of continuous point processes (vecplay)
  // e.g. stimulus. vecplay->pd points to a read-only var used by
  // point proc mechanisms' nrn_current function
  nt->_t += nt->_dt * 0.5;  // Backward-Euler half-step
  branch->FixedPlayContinuous();

  // Sets RHS an D to zero
  solver::HinesSolver::ResetMatrixRHSandD(branch);

  // current sums i and didv to parent ion, and adds contribnutions to RHS and D
  branch->CallModFunction(Mechanism::ModFunctions::kCurrent);

  // add parent and children currents (A*dv and B*dv) to RHS
  solver::HinesSolver::SetupMatrixRHS(branch);

  // update positions holding jacobians (does nothing so far)
  branch->CallModFunction(Mechanism::ModFunctions::kJacob);
  //(so far only calls nrn_jacob_capacitance, which sums contributions to D)
  branch->CallModFunction(Mechanism::ModFunctions::kJacobCapacitance);

  // backup up D (will be used for jacobian diagonal)
  memcpy(user_data->jacob_d_, nt->_actual_d, sizeof(double) * nt->end);

  // add parent and children currents (A and B) to D
  solver::HinesSolver::SetupMatrixDiagonal(branch);

  // Gaussian Elimination (sets dV/dt=RHS[i])
  branch->SolveTreeMatrix();

  // update ions currents based on RHS and dI/dV
  second_order_cur(branch->nt_, input_params_->second_order_);

  // update capacitance current based on RHS and dI/dV
  branch->CallModFunction(Mechanism::ModFunctions::kCurrentCapacitance);

  // updates V: v[i] += second_order_multiplier * rhs[i]
  solver::HinesSolver::UpdateV(branch);

  // update mechanisms state (eg opening vars) based on voltage
  // State is now given by CVODE
  // branch->CallModFunction(Mechanism::ModFunctions::kState);

  // compute ydot=f(t,y) for values of mechanisms opening vars
  branch->CallModFunction(Mechanism::ModFunctions::kODEMatsol);

  // set ydot with RHS state, and recover NrnThread->data state
  CvodesAlgorithm::CopyRHSToYdot(branch, ydot);

  nt->_t = t;
  return CV_SUCCESS;
}

// jacobian routine: compute J(t,y) = df/dy
int CvodesAlgorithm::JacobianFunction(long int N, realtype t, N_Vector y,
                                      N_Vector fy, DlsMat J,
                                      void *user_data_ptr, N_Vector, N_Vector,
                                      N_Vector) {
  realtype **jac = J->cols;
  BranchCvodes::UserData *user_data = (BranchCvodes::UserData *)user_data_ptr;
  Branch *branch = user_data->branch_;
  const NrnThread *nt = branch->nt_;
  assert(t == user_data->rhs_last_time_);

  const int a_offset = nt->_actual_a - nt->_data;
  const int b_offset = nt->_actual_b - nt->_data;
  const int compartments_count = nt->end;

  const double *a = &nt->_data[a_offset];
  const double *b = &nt->_data[b_offset];
  const double *d = user_data->jacob_d_;  // computed by RHS
  const int *p = nt->_v_parent_index;

  // Jacobian for main current equation:
  // d(dV/dt) / dV     = sum_i g_i x_i  + D  //mechs currents + D
  // d(dV/dt) / dV_p   = -A          //parent compartment
  // d(dV/dt) / dV_c_i = -B_c_i  //children_i compartment

  // jac[a][b] = d/dV_b (dV_a/dt)
  for (int n = 0; n < compartments_count; n++) {
    // if not stepping backwards
    // assert(user_data->rhs_last_time_ != user_data->rhs_second_last_time_);
    if (user_data->rhs_last_time_ > user_data->rhs_second_last_time_) {
      assert(d[n] >= 0);  // positive (currents and mechs contribution, if any)
      assert(a[n] <= 0 && b[n] <= 0);  // negative (resistance)
    }
    jac[n][n] = d[n];  // C = d/dV_n (dV_n/dt)
    if (n == 0) continue;
    jac[p[n]][n] = b[n];  // B = d/dV_n (dV_p/dt)
    jac[n][p[n]] = a[n];  // A = d/dV_p (dV_n/dt)
  }

  // get new derivative of mechs states (ions do not have,
  // capacitance added already to D in nrn_jacob_capacitance)
  branch->CallModFunction(Mechanism::ModFunctions::kODESpec);

  int state_dv_index=-1;
  int compartment_id=-1;
  int jac_offset = compartments_count;
  for (int m = 0; m < neurox::mechanisms_count_; m++) {
    Mechanism *mech = mechanisms_[m];
    Memb_list *mech_instances = &branch->mechs_instances_[m];
    double * mech_instances_data = mech_instances->data;
    int ml_data_offset = 0;
    for (int n = 0; n < mech_instances->nodecount; n++) {
      for (int s = 0; s < mech->state_vars_->count_; s++) {
        state_dv_index = mech->state_vars_->dv_offsets_[s];
        compartment_id = mech_instances->nodeindices[n];
#if LAYOUT == 1
        int state_dv_offset =
            ml_data_offset + mech->data_size_ * n + state_var_dv_index;
#else
        int state_dv_offset =
            ml_data_offset +
            tools::Vectorizer::SizeOf(mech_instances->nodecount) * state_dv_index +
            n;
#endif
        jac[jac_offset][compartment_id] = mech_instances_data[state_dv_offset];
        jac_offset++;
      }
      ml_data_offset += tools::Vectorizer::SizeOf(mech_instances->nodecount) *
                        mech->data_size_;
    }
  }

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
      user_data_(nullptr),
      y_(nullptr),
      spikes_lco_(HPX_NULL) {}

CvodesAlgorithm::BranchCvodes::~BranchCvodes() {
  N_VDestroy_Serial(y_);   /* Free y vector */
  CVodeFree(&cvodes_mem_); /* Free integrator memory */
  delete this->user_data_;
}

CvodesAlgorithm::BranchCvodes::UserData::UserData(Branch *branch)
    : branch_(branch), rhs_last_time_(0.0) {
  this->jacob_d_ = new double[branch->nt_->end];
  this->data_bak_ = new double[branch->nt_->_ndata];
}

CvodesAlgorithm::BranchCvodes::UserData::~UserData() {
  delete[] data_bak_;
  delete[] jacob_d_;
  branch_ = nullptr;
}

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
  int map_offset = 0;
  branch_cvodes->state_var_map_ =
      new double *[equations_count - compartments_count];
  branch_cvodes->state_dv_map_ =
      new double *[equations_count - compartments_count];
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
            tools::Vectorizer::SizeOf(mech_instances->nodecount) * state_var_index +
            n;

        int state_dv_offset =
            ml_data_offset +
            tools::Vectorizer::SizeOf(mech_instances->nodecount) * state_dv_index +
            n;
#endif
        branch_cvodes->state_var_map_[map_offset] =
            &mech_instances->data[state_var_offset];
        branch_cvodes->state_dv_map_[map_offset] =
            &mech_instances->data[state_dv_offset];
        map_offset++;
      }
    }
    ml_data_offset += tools::Vectorizer::SizeOf(mech_instances->nodecount) *
                        mech->data_size_;
  }
  branch_cvodes->y_ = N_VMake_Serial(equations_count, y_data);

  // absolute tolerance array (low for voltages, high for mech states)
  branch_cvodes->absolute_tolerance_ = N_VNew_Serial(equations_count);
  for (int i = 0; i < compartments_count; i++)
    NV_Ith_S(branch_cvodes->absolute_tolerance_, i) = kAbsToleranceVoltage;
  for (int i = compartments_count; i < equations_count; i++)
    NV_Ith_S(branch_cvodes->absolute_tolerance_, i) = kAbsToleranceMechStates;

  // CVodeCreate creates an internal memory block for a problem to
  // be solved by CVODES, with Backward Differentiation (or Adams)
  // and Newton solver (recommended for stiff problems, see header)
  cvodes_mem = CVodeCreate(CV_BDF, CV_NEWTON);

  // CVodeInit allocates and initializes memory for a problem
  flag = CVodeInit(cvodes_mem, CvodesAlgorithm::RHSFunction, 0.0 /*t_0*/,
                   branch_cvodes->y_);
  assert(flag == CV_SUCCESS);

  // specify integration tolerances. MUST be called before CVode.
  flag = CVodeSVtolerances(cvodes_mem, kRelativeTolerance,
                           branch_cvodes->absolute_tolerance_);
  assert(flag == CV_SUCCESS);

  // specify user data to be used on functions as void* user_data_ptr;
  branch_cvodes->user_data_ = new BranchCvodes::UserData(local);
  flag = CVodeSetUserData(cvodes_mem, branch_cvodes->user_data_);
  assert(flag == CV_SUCCESS);

  // specify g as root function and roots
  int roots_direction[1] = {1};  // AP threshold reached for increasing voltage
  flag = CVodeRootInit(cvodes_mem, 1, CvodesAlgorithm::RootFunction);
  assert(flag == CV_SUCCESS);

// initializes the memory record and sets various function
// fields specific to the dense linear solver module.
// Note: direct solvers give the solution (LU-decomposition, etc)
// Indirect solvers require iterations (eg Jacobi method)

#if NEUROX_CVODES_JACOBIAN_SOLVER == 0  // dense solver
  flag = CVDense(cvodes_mem, equations_count);
  flag = CVDlsSetDenseJacFn(cvodes_mem, CvodesAlgorithm::JacobianFunction);
#else  // sparse colver
  // Requires installation of Superlumt or KLU
  flag =
      CVSlsSetSparseJacFn(cvode_mem_, CvodesAlgorithm::JacobianSparseFunction);
  int nnz = equations_count * equations_count;
#if NEUROX_CVODES_JACOBIAN_SOLVER == 1
  flag = CVKLU(cvode_mem, 1, equations_count, nnz);
#else  //==2
  flag = CVSuperLUMT(cvode_mem, 1, equations_count, nnz);
#endif
  assert(flag == CV_SUCCESS);
#endif

  assert(flag == CV_SUCCESS);

  // specify the dense (user-supplied) Jacobian function. Compute J(t,y).
  // see chapter 8 -- providing alternate linear solver modules
  assert(flag == CV_SUCCESS);

  // TODO
  // CVodeSetInitStep(cvodes_mem, kMinStepSize);
  // CVodeSetMinStep(cvodes_mem, kMinStepSize);
  CVodeSetMaxStep(cvodes_mem,
                  CoreneuronAlgorithm::CommunicationBarrier::kCommStepSize);
  CVodeSetStopTime(cvodes_mem, input_params_->tstop_);
  CVodeSetMaxOrd(cvodes_mem, 5);  // max order of the BDF method
  CVodeSetRootDirection(cvodes_mem, roots_direction);

  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t CvodesAlgorithm::BranchCvodes::Run = 0;
int CvodesAlgorithm::BranchCvodes::Run_handler() {
  NEUROX_MEM_PIN(neurox::Branch);
  assert(local->soma_);
  BranchCvodes *branch_cvodes =
      (BranchCvodes *)local->soma_->algorithm_metadata_;
  NrnThread *nt = local->nt_;
  void *cvodes_mem = branch_cvodes->cvodes_mem_;

  int roots_found[1];  // AP threshold
  int flag = CV_ERR_FAILURE;
  realtype tout = 0;

  while (nt->_t < input_params_->tstop_) {
    // delivers all events whithin the next min step size
    local->DeliverEvents(nt->_t + CvodesAlgorithm::kEventsDeliveryTimeWindow);

    // get tout as time of next undelivered event (if any)
    hpx_lco_sema_p(local->events_queue_mutex_);
    if (!local->events_queue_.empty()) tout = local->events_queue_.top().first;
    hpx_lco_sema_v_sync(local->events_queue_mutex_);
    tout = std::min(input_params_->tstop_, tout);

    // call CVODE method: steps until reaching/passing tout;
    // TODO can it walk backwards for missed event?
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
