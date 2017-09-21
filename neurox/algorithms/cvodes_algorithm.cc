#include "neurox/algorithms/cvodes_algorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

double CvodesAlgorithm::BranchCvodes::min_step_size_ = -1;

CvodesAlgorithm::CvodesAlgorithm(){}

CvodesAlgorithm::~CvodesAlgorithm() {}

const AlgorithmId CvodesAlgorithm::GetId() {
  return AlgorithmId::kCvodes;
}

const char* CvodesAlgorithm::GetString() {
  return "CVODES";
}

void CvodesAlgorithm::UpdateNrnThreadFromCvodeState(Branch * branch, N_Vector y)
{
    BranchCvodes* branch_cvodes = (BranchCvodes*) branch->soma_->algorithm_metadata_;
    double * y_data = NV_DATA_S(y);

    //copy voltage
    int compartments_count = branch->nt_->end;
    memcpy(branch->nt_->_actual_v, y_data, sizeof(double)*compartments_count);

    //use map to copy states to NrnThread->data
    for (int i=0; i<branch_cvodes->equations_count_ - compartments_count; i++)
        *branch_cvodes->equations_map_[i] = y_data[i];
}

/// g root function to compute g_i(t,y) .
int CvodesAlgorithm::RootFunction(realtype t, N_Vector y, realtype *gout, void *user_data_ptr)
{
    Branch * branch = ((BranchCvodes::UserData*) user_data_ptr)->branch_;
    int v_offset = branch->thvar_ptr_- branch->nt_->_actual_v;
    double v = NV_Ith_S(y, v_offset);

    //How it works: when gout[x] is zero, a root is found
    assert(v>=-100 && v<30);
    gout[0] = v - branch->soma_->threshold_; //AP threshold reached
    return CV_SUCCESS;
}

int CvodesAlgorithm::JacobianSparseMatrix(realtype t,
               N_Vector y, N_Vector fy, SlsMat JacMat, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    assert(0);
    return(0);
}

/// f routine to compute ydot=f(t,y), i.e. compute new values of nt->data
int CvodesAlgorithm::RHSFunction(realtype t, N_Vector y, N_Vector ydot, void *user_data_ptr)
{
    BranchCvodes::UserData *user_data = (BranchCvodes::UserData*)  user_data_ptr;
    Branch * branch = user_data->branch_;
    NrnThread * nt = branch->nt_;

    //NOTE: status has to be full recoverable as it will step with several ts,
    //and update the status based on the best f(y,t) found;
    //Therefore, we back up NrnThread->data and time (weights are only
    //changed on net_receive), and we recoved them later

    //backup previous state
    double t_bak = nt->_t;
    memcpy(user_data->data_bak_, nt->_data, sizeof(double)*nt->_ndata);

    //Note: due to current stae of the art, this function will compute
    //both the RHS  (vec_RHS) and jacobian (D).

    //update vars in NrnThread->data described by our CVODES state
    CvodesAlgorithm::UpdateNrnThreadFromCvodeState(branch, y);

    /////////   Get new RHS and D from current state ///////////

    nt->_t += .5 * BranchCvodes::min_step_size_;

    //Updates internal states of continuous point processes (vecplay)
    //e.g. stimulus. vecplay->pd points to a read-only var used by
    //point proc mechanisms' nrn_current function
    branch->FixedPlayContinuous();

    // Sets RHS an D to zero
    solver::HinesSolver::ResetMatrixRHSandD(branch);

    // current sums i and didv to parent ion, and adds contribnutions to RHS and D
    branch->CallModFunction(Mechanism::ModFunctions::kCurrent);

    //add parent and children currents (A*dv and B*dv) to RHS
    solver::HinesSolver::SetupMatrixRHS(branch);

    //update positions holding jacobians (does nothing so far)
    branch->CallModFunction(Mechanism::ModFunctions::kJacob);
    //(so far only calls nrn_jacob_capacitance, which sums contributions to D)
    branch->CallModFunction(Mechanism::ModFunctions::kJacobCapacitance);

    //backup up D (will be used for jacobian diagonal)
    memcpy(user_data->jacob_d_, nt->_actual_d, sizeof(double)*nt->end);

    //add parent and children currents (A and B) to D
    solver::HinesSolver::SetupMatrixDiagonal(branch);

    //Gaussian Elimination (sets dV/dt=RHS[i])
    branch->SolveTreeMatrix();

    //update ions currents based on RHS and dI/dV
    second_order_cur(branch->nt_, input_params_->second_order_);

    //update capacitance current based on RHS and dI/dV
    branch->CallModFunction(Mechanism::ModFunctions::kCurrentCapacitance);

    //updates V: v[i] += second_order_multiplier * rhs[i]
    solver::HinesSolver::UpdateV(branch);

    nt->_t += .5 * BranchCvodes::min_step_size_;

    // update mechanisms state (eg opening vars) based on voltage
    //TODO until we add mechanisms, changes below will be ignored, as state is overwritten
    branch->CallModFunction(Mechanism::ModFunctions::kState);

    //set ydot with RHS state
    memcpy(NV_DATA_S(ydot), branch->nt_->_actual_rhs, sizeof(branch->nt_->end));

    //recover previous state
    nt->_t = t_bak;
    memcpy(nt->_data, user_data->data_bak_, sizeof(double)*nt->_ndata);

    return CV_SUCCESS;
}

//jacobian routine: compute J(t,y) = df/dy
int CvodesAlgorithm::JacobianFunction(
        long int N, realtype t,
        N_Vector y, N_Vector fy,
        DlsMat J, void *user_data_ptr,
        N_Vector, N_Vector, N_Vector)
{
    realtype ** jac = J->cols;
    BranchCvodes::UserData *user_data = (BranchCvodes::UserData*)  user_data_ptr;
    const Branch * branch = user_data->branch_;
    const BranchCvodes* branch_cvodes = (BranchCvodes*) branch->soma_->algorithm_metadata_;
    assert(N==branch_cvodes->equations_count_);

    const int a_offset = branch->nt_->_actual_a - branch->nt_->_data;
    const int b_offset = branch->nt_->_actual_b - branch->nt_->_data;
    const int compartments_count = branch->nt_->end;

    const double *a = &branch->nt_->_data[a_offset];
    const double *b = &branch->nt_->_data[b_offset];
    const double *d = user_data->jacob_d_;
    const int *p = branch->nt_->_v_parent_index;

    //Jacobian for main current equation:
    // d(C dV/dt) / dV     = sum_i g_i x_i  + D  //mechs currents + D
    // d(C dV/dt) / dV_p   = -A          //parent compartment
    // d(C dV/dt) / dV_c_i = -B_c_i  //children_i compartment

    // simillar to HinesSolver::SetupMatrixDiagonal
    // Reminder: vec_d is the derivative of V so we sum
    // partial derivatives A and B for parents/children contribution
    for (int i=0; i<compartments_count; i++)
    {
        jac[i][i]    = d[i]; //diagonal
        if (i==0) continue;
        jac[i][p[i]] = b[i]; //lower-diag (children) in same row
        jac[p[i]][i] = a[i]; //upper-diag (parents) in same column
    }

    /*
    //USE voltage to update states
    int v_offset = branch->nt_->_actual_v - branch->nt_->_data;
    realtype *v = NV_DATA_S(y);

    //Add di/dv contributions from mechanisms currents to current equation
    int compartment_id=-1, y_data_offset=-1, g_data_index=-1;
    int data_offset=compartments_count;
    realtype y_val=-1;
    const double cfac = .001 * branch->nt_->cj; //capac.c::nrn_jacob_capacitance
    for (int m = 0; m < neurox::mechanisms_count_; m++)
    {
        Mechanism * mech=mechanisms_[m];
        Memb_list* mech_instances = &branch->mechs_instances_[m];

        if (mech->memb_func_.current != NULL) //if it has di/dv
        {
          for (int n=0; n<mech_instances->nodecount; n++)
          {
            g_data_index = mech->state_vars_offsets_[0];
#if LAYOUT == 1
            y_data_offset = data_offset + mech->data_size_ * n + g_data_index;
#else
            y_data_offset = data_offset + tools::Vectorizer::SizeOf(mech_instances->nodecount) * g_data_index + n;
#endif
            y_val = v[y_data_offset];
            if (mech->is_ion_) //y is second order current (mechs functions di/dV)
            {
                //TODO this is wrong
              //cur variable is one position before dcurdv (eion.c)
              jac[y_data_offset-1][v_offset+compartment_id] = y_val; //dcurdv
            }
            else //y is first order current (current function vec_D for currents C*dV/dt)
            {
              //NOTE: vec_D is d(dV/dt) so y_val is the jacob of dV/dt
              if (mech->type_==MechanismTypes::kCapacitance)
                  y_val*=cfac; //eion.c::nrn_jacob_capacitance
              compartment_id = mech_instances->nodeindices[m];
              assert(jac[compartment_id][y_data_offset]==0);
              jac[compartment_id][y_data_offset] = y_val;


              //TODO aren't we missing shadow_didv
              //dcurdv for current mechanism
            }
          }
        }
        data_offset += tools::Vectorizer::SizeOf(mech_instances->nodecount)*mech->data_size_;
    }
    */
    return CV_SUCCESS;
}

/////////////////////// Algorithm abstract class ////////////////////////


void CvodesAlgorithm::Init() {
  BranchCvodes::min_step_size_ = input_params_->dt_ * 0.1;
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


void CvodesAlgorithm::StepBegin(Branch*){
}

void CvodesAlgorithm::StepEnd(Branch* b, hpx_t spikesLco) {
}

void CvodesAlgorithm::Run(Branch* b, const void* args){
}

hpx_t CvodesAlgorithm::SendSpikes(Neuron* n, double tt, double) {
    return Neuron::SendSpikesAsync(n, tt);
}


//////////////////////////// BranchCvodes /////////////////////////

CvodesAlgorithm::BranchCvodes::BranchCvodes()
    :cvodes_mem_(nullptr), equations_count_(-1),
     equations_map_(nullptr),  user_data_(nullptr),
     y_(nullptr)
{}

CvodesAlgorithm::BranchCvodes::~BranchCvodes()
{
    N_VDestroy_Serial(y_); /* Free y vector */
    CVodeFree(&cvodes_mem_); /* Free integrator memory */
    delete this->user_data_;
}

CvodesAlgorithm::BranchCvodes::UserData::UserData(Branch * branch)
{
    this->branch_ = branch;
    this->jacob_d_ = new double [branch->nt_->end];
    this->data_bak_ = new double [branch->nt_->_ndata];
}

CvodesAlgorithm::BranchCvodes::UserData::~UserData()
{
    delete [] data_bak_;
    delete [] jacob_d_;
    branch_=nullptr;
}

hpx_action_t CvodesAlgorithm::BranchCvodes::Init = 0;
int CvodesAlgorithm::BranchCvodes::Init_handler()
{
    NEUROX_MEM_PIN(neurox::Branch);
    assert(local->soma_->algorithm_metadata_);
    BranchCvodes * branch_cvodes = (BranchCvodes*) local->soma_->algorithm_metadata_;
    void *& cvodes_mem = branch_cvodes->cvodes_mem_;
    int & equations_count = branch_cvodes->equations_count_;
    int compartments_count = local->nt_->end;
    NrnThread *& nt = local->nt_;

    //set dt to 1, so that dV/dt equations are on the right scale
    nt->_dt=1.0;

    int flag = CV_ERR_FAILURE;

    //equations: voltages per compartments + mechanisms * states
    equations_count = compartments_count;
    /*
    for (int m = 0; m < neurox::mechanisms_count_; m++)
    {
        Mechanism * mech=mechanisms_[m];
        Memb_list* mech_instances = &local->mechs_instances_[m];
        equations_count += mech_instances->nodecount * mech->state_vars_count_;
    }
    */

    //create initial state y_0 for state array y
    floble_t * y_data = new floble_t[equations_count];
    for (int i=0; i<compartments_count; i++)
      y_data[i] = nt->_actual_v[i];

    //create map from y to NrnThread->data (mech-states)
    branch_cvodes->equations_map_ = new double*[equations_count - compartments_count];
    int equations_map_offset=0;
    for (int m = 0; m < neurox::mechanisms_count_; m++)
    {
        Mechanism * mech = mechanisms_[m];
        Memb_list* mech_instances = &local->mechs_instances_[m];
        int ml_data_offset=0;
        for (int n=0; n<mech_instances->nodecount; n++)
        {
            for (int s=0; s<mech->state_vars_count_; s++)
            {
              int state_index = mech->state_vars_offsets_[s];
#if LAYOUT == 1
              int state_data_offset = ml_data_offset + mech->data_size_ * n +  state_index;
#else
              int state_data_offset = ml_data_offset + tools::Vectorizer::SizeOf(mech_instances->nodecount) * state_index + n;
#endif
              branch_cvodes->equations_map_[equations_map_offset] = &mech_instances->data[state_data_offset];
              equations_map_offset++;
            }
            ml_data_offset += tools::Vectorizer::SizeOf(mech_instances->nodecount) * mech->data_size_;
        }
        equations_map_offset += mech_instances->nodecount * mech->state_vars_count_;
    }
    branch_cvodes->y_ = N_VMake_Serial(equations_count, y_data);

    //absolute tolerance array (low for voltages, high for mech states)
    branch_cvodes->absolute_tolerance_ = N_VNew_Serial(equations_count);
    for (int i=0; i<compartments_count; i++)
        NV_Ith_S(branch_cvodes->absolute_tolerance_, i) = kAbsToleranceVoltage;
    for (int i=compartments_count; i<equations_count; i++)
        NV_Ith_S(branch_cvodes->absolute_tolerance_, i) = kAbsToleranceMechStates;

    //CVodeCreate creates an internal memory block for a problem to
    //be solved by CVODES, with Backward Differentiation (or Adams)
    //and Newton solver (recommended for stiff problems, see header)
    cvodes_mem = CVodeCreate(CV_BDF, CV_NEWTON);

    //CVodeInit allocates and initializes memory for a problem
    flag = CVodeInit(cvodes_mem, CvodesAlgorithm::RHSFunction, 0.0 /*t_0*/, branch_cvodes->y_);
    assert(flag==CV_SUCCESS);

    //specify integration tolerances. MUST be called before CVode.
    flag = CVodeSVtolerances(cvodes_mem, kRelativeTolerance, branch_cvodes->absolute_tolerance_);
    assert(flag==CV_SUCCESS);

    //specify user data to be used on functions as void* user_data_ptr;
    branch_cvodes->user_data_ = new BranchCvodes::UserData(local);
    flag = CVodeSetUserData(cvodes_mem,  branch_cvodes->user_data_);
    assert(flag==CV_SUCCESS);

    //specify g as root function and roots
    int roots_direction [1] = {1}; //AP threshold reached for increasing voltage
    flag = CVodeRootInit(cvodes_mem, 1, CvodesAlgorithm::RootFunction);
    assert(flag==CV_SUCCESS);

    //initializes the memory record and sets various function
    //fields specific to the dense linear solver module.
    //Note: direct solvers give the solution (LU-decomposition, etc)
    //Indirect solvers require iterations (eg Jacobi method)

#if NEUROX_CVODES_JACOBIAN_SOLVER==0 //dense solver
    flag = CVDense(cvodes_mem, equations_count);
    flag = CVDlsSetDenseJacFn(cvodes_mem, CvodesAlgorithm::JacobianFunction);
#else //sparse colver
    //Requires installation of Superlumt or KLU
    flag = CVSlsSetSparseJacFn(cvode_mem_, CvodesAlgorithm::JacobianSparseFunction);
    int nnz = equations_count * equations_count;
#if NEUROX_CVODES_JACOBIAN_SOLVER==1
    flag = CVKLU(cvode_mem, 1, equations_count, nnz);
#else //==2
    flag = CVSuperLUMT(cvode_mem, 1, equations_count, nnz);
#endif
    assert(flag==CV_SUCCESS);
#endif

    assert(flag==CV_SUCCESS);

    //specify the dense (user-supplied) Jacobian function. Compute J(t,y).
    //see chapter 8 -- providing alternate linear solver modules
    assert(flag==CV_SUCCESS);

    CVodeSetInitStep(cvodes_mem, min_step_size_);
    CVodeSetMinStep(cvodes_mem, min_step_size_);
    CVodeSetMaxStep(cvodes_mem, CoreneuronAlgorithm::CommunicationBarrier::kCommStepSize);
    CVodeSetStopTime(cvodes_mem, input_params_->tstop_);
    CVodeSetMaxOrd(cvodes_mem, 5); //max order of the BDF method
    CVodeSetRootDirection(cvodes_mem, roots_direction);

    return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t CvodesAlgorithm::BranchCvodes::Run = 0;
int CvodesAlgorithm::BranchCvodes::Run_handler()
{
    NEUROX_MEM_PIN(neurox::Branch);
    assert(local->soma_);
    BranchCvodes * branch_cvodes = (BranchCvodes*) local->soma_->algorithm_metadata_;
    NrnThread * nt = local->nt_;
    void * cvodes_mem = branch_cvodes->cvodes_mem_;

    int roots_found[1]; //AP threshold

    int flag = CV_ERR_FAILURE;
    realtype tout = 0;
    hpx_t spikes_lco = HPX_NULL;

    while(nt->_t < input_params_->tstop_)
    {
      //get tout as time of next undelivered event (if any)
      hpx_lco_sema_p(local->events_queue_mutex_);
      if (!local->events_queue_.empty())
          tout = local->events_queue_.top().first;
      hpx_lco_sema_v_sync(local->events_queue_mutex_);
      tout = std::min(input_params_->tstop_, tout);

      //call CVODE method: steps until reaching/passing tout;
      //TODO can it walk backwards for missed event?
      flag = CVode(cvodes_mem, tout, branch_cvodes->y_, &(nt->_t), CV_NORMAL);

      printf("At t = %0.4e   V =%14.6e  %14.6e  %14.6e\n",
             nt->_t, NV_Ith_S(branch_cvodes->y_,0),
             NV_Ith_S(branch_cvodes->y_,1), NV_Ith_S(branch_cvodes->y_,2));

      if(flag==CV_ROOT_RETURN) //CVODE succeeded and roots found
      {
        //delivers all events whithin the next min step size
        local->DeliverEvents(nt->_t + BranchCvodes::min_step_size_);

        flag = CVodeGetRootInfo(cvodes_mem, roots_found);
        assert(flag==CV_SUCCESS);
        assert(roots_found[0]!=0); //only possible root: AP threshold reached
        //NOTE: if root is found, integrator time is now at time of root!
        //(+1 value ascending, -1 valued descending)
        if (roots_found[0] > 0) //AP threshold reached from below
        {
           //TODO we can use a root to stop progress in time and wait for deliver?
           //use several hpx-all reduce to mark performance?
           spikes_lco=local->soma_->SendSpikes(nt->_t);
        }
      }

      //Success from CVodeGetRootInfo or CVode
      if (flag == CV_SUCCESS)
        nt->_t = tout;

      /*
      long num_steps=-1, num_jacob_evals=-1, num_rhs_evals=-1;
#if NEUROX_CVODES_JACOBIAN_SOLVER==0
      CVDlsGetNumJacEvals(cvodes_mem, &num_jacob_evals);
      CVDlsGetNumRhsEvals(cvodes_mem, &num_rhs_evals);
#else
      CVSlsGetNumJacEvals(cvodes_mem, &num_jacob_evals);
      CVSlsGetNumRhsEvals(cvodes_mem, &num_rhs_evals);
#endif

      CVodeGetNumSteps(cvodes_mem, &num_steps);
      realtype last_step_size=-1;
      CVodeGetLastStep(cvodes_mem, &last_step_size);
      N_Vector estimated_local_errors = N_VNewEmpty_Serial(branch_cvodes->equations_count_);
      CVodeGetEstLocalErrors(cvodes_mem, estimated_local_errors);
      */
    }

    /* Print some final statistics */
    //PrintFinalStats(cvode_mem_);

    return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t CvodesAlgorithm::BranchCvodes::Clear = 0;
int CvodesAlgorithm::BranchCvodes::Clear_handler()
{
    NEUROX_MEM_PIN(neurox::Branch);
    assert(local->soma_);
    BranchCvodes * branch_cvodes = (BranchCvodes*) local->soma_->algorithm_metadata_;
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
