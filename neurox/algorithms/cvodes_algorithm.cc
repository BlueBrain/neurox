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

/// f routine to compute ydot=f(t,y), i.e. compute new values of nt->data
void CvodesAlgorithm::ReevaluateBranch(Branch * branch)
{
    //delivers all events whithin the next min step size
    branch->DeliverEvents(branch->nt_->_t + BranchCvodes::min_step_size_);

    //Updates internal states of continuous point processes (vecplay)
    //e.g. stimulus. vecplay->pd points to a read-only var used by
    //point proc mechanisms' nrn_current function
    branch->FixedPlayContinuous();

    ////// From HinesSolver::SetupTreeMatrix ///////

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

    //add parent and children currents (A and B) to D
    solver::HinesSolver::SetupMatrixDiagonal(branch);


    ////// From HinesSolver::SolveTreeMatrix ///////
    //Gaussian Elimination (sets RHS)
    branch->SolveTreeMatrix();


    ////// From eion.c::second_order_cur ///////
    //updates ionic currents as rhs*didv from the children
    second_order_cur(branch->nt_, input_params_->second_order_);


    ////// From HinesSolver::UpdateV ///////
    // v[i] += second_order_multiplier * rhs[i];
    solver::HinesSolver::UpdateV(branch);

    branch->CallModFunction(Mechanism::ModFunctions::kCurrentCapacitance);

    ////// From main loop - call state function
    branch->CallModFunction(Mechanism::ModFunctions::kState);
}


void CvodesAlgorithm::UpdateNrnThreadFromCvodeState(Branch * branch)
{
    BranchCvodes* branch_cvodes = (BranchCvodes*) branch->soma_->algorithm_metadata_;
    double * y = NV_DATA_S(branch_cvodes->y_);
    for (int i=0; i<branch_cvodes->equations_count_; i++)
        *branch_cvodes->equations_map_[i] = y[i];
}

/// g root function to compute g_i(t,y) .
int CvodesAlgorithm::RootFunction(realtype t, N_Vector y, realtype *gout, void *user_data)
{
    Branch * b = (Branch*) user_data;
    realtype v = *b->thvar_ptr_;

    //How it works: when gout[x] is zero, a root is found
    gout[0] = v - b->soma_->threshold_; //AP threshold reached
#ifndef NDEBUG
    gout[1] = v - 30; //Debug: reached   50 mV (too high)
    gout[2] = v + 90; //Debug: reached -100 mV (too low )
#endif
    return 0;
}

int CvodesAlgorithm::JacobianSparseMatrix(realtype t,
               N_Vector y, N_Vector fy, SlsMat JacMat, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    assert(0);
    return(0);
}

int CvodesAlgorithm::RHSFunction(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    Branch * branch = (Branch*) user_data;

    int a_offset = branch->nt_->_actual_a - branch->nt_->_data;
    int b_offset = branch->nt_->_actual_b - branch->nt_->_data;
    const double *a = &branch->nt_->_data[a_offset];
    const double *b = &branch->nt_->_data[b_offset];
    int * p = branch->nt_->_v_parent_index;
    int compartments_count = branch->nt_->end;

    double * v   = NV_DATA_S(y);
    double * vdot = NV_DATA_S(ydot);

    // contribution from parent and children compartments
    // d/dV (C dV/dt) = b * V_p + sum_i a_i V_c_i (parent + sum of all children)
    // reminder: a and b arrays are constant
    // Note: from "Coreneuron Overview" report from Ben (page 4)
    // vdot[i] is rhs[i] which is the r.h.s. of dV/dt= C * ...

    for (int i = 0; i < compartments_count; i++) {
      assert(vdot[i]==0); //checking if its necessary or not
      vdot[i] = 0;
    }

    // simillar to HinesSolver::SetupMatrixRHS
    double dv=-1;
    for (int i=0; i<compartments_count; i++)
    {
        dv = v[p[i]] - v[i];  // reads from parent
        vdot[i] -= b[i] * dv;
        vdot[p[i]] += a[i] * dv;  // writes to parent    }
    }
}

//jacobian routine: compute J(t,y) = df/dy
int CvodesAlgorithm::JacobianFunction(
        long int N, realtype t,
        N_Vector y, N_Vector fy,
        DlsMat J, void *user_data,
        N_Vector, N_Vector, N_Vector)
{
    realtype ** jac = J->cols;
    Branch * branch = (Branch*) user_data;
    BranchCvodes* branch_cvodes = (BranchCvodes*) branch->soma_->algorithm_metadata_;
    assert(N==branch_cvodes->equations_count_);

    int a_offset = branch->nt_->_actual_a - branch->nt_->_data;
    int b_offset = branch->nt_->_actual_b - branch->nt_->_data;
    int v_offset = branch->nt_->_actual_v - branch->nt_->_data;
    int compartments_count = branch->nt_->end;

    double *a = &branch->nt_->_data[a_offset];
    double *b = &branch->nt_->_data[b_offset];
    realtype *v = NV_DATA_S(y);
    int *p = branch->nt_->_v_parent_index;

    //Jacobian for main current equation:
    // d(C dV/dt) / dV     = sum_i g_i x_i   //mechs currents
    // d(C dV/dt) / dV_p   = -A          //parent compartment
    // d(C dV/dt) / dV_c_i = -B_c_i  //children_i compartment

    // simillar to HinesSolver::SetupMatrixDiagonal
    // Reminder: vec_d is the derivative of V so we sum
    // partial derivatives A and B for parents/children contribution
    for (int i=0; i<compartments_count; i++)
    {
        jac[i][p[i]] -= a[i];
        jac[p[i]][i] -= b[i];
    }

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

    return 0;
}

/////////////////////// Algorithm abstract class ////////////////////////


void CvodesAlgorithm::Init() {
  BranchCvodes::min_step_size_ = input_params_->dt_ * 0.5;
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
    :cvodes_mem_(nullptr), equations_count_(-1)
{}

CvodesAlgorithm::BranchCvodes::~BranchCvodes()
{
    N_VDestroy_Serial(y_); /* Free y vector */
    CVodeFree(&cvodes_mem_); /* Free integrator memory */
}

hpx_action_t CvodesAlgorithm::BranchCvodes::Init = 0;
int CvodesAlgorithm::BranchCvodes::Init_handler()
{
    NEUROX_MEM_PIN(neurox::Branch);
    assert(local->soma_);
    BranchCvodes * branch_cvodes = (BranchCvodes*) local->soma_->algorithm_metadata_;
    N_Vector absolute_tolerance = branch_cvodes->absolute_tolerance_;
    void * cvodes_mem = branch_cvodes->cvodes_mem_;
    int compartments_count = local->nt_->end;
    NrnThread *& nt = local->nt_;

    int flag = CV_ERR_FAILURE;

    //equations: voltages per compartments + mechanisms * states
    int & equations_count = branch_cvodes->equations_count_ = 0;
    equations_count += compartments_count;
    for (int m = 0; m < neurox::mechanisms_count_; m++)
    {
        Mechanism * mech=mechanisms_[m];
        Memb_list* mech_instances = &local->mechs_instances_[m];
        equations_count += mech_instances->nodecount * mech->state_vars_count_;
    }

    //create initial state y_0 for state array y
    floble_t * y_data = new floble_t[equations_count];
    for (int i=0; i<compartments_count; i++)
    {
      y_data[i] = nt->_actual_v[i];
      branch_cvodes->equations_map_[i]=&local->nt_->_actual_v[i];
    }

    //create map from y to NrnThread->data (mech-states)
    branch_cvodes->equations_map_ = new double*[equations_count];
    int equations_map_offset=compartments_count;
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
    absolute_tolerance = N_VNew_Serial(equations_count);
    for (int i=0; i<compartments_count; i++)
        NV_Ith_S(absolute_tolerance, i) = kAbsToleranceVoltage;
    for (int i=compartments_count; i<equations_count; i++)
        NV_Ith_S(absolute_tolerance, i) = kAbsToleranceMechStates;

    //CVodeCreate creates an internal memory block for a problem to
    //be solved by CVODES, with Backward Differentiation (or Adams)
    //and Newton solver (recommended for stiff problems, see header)
    cvodes_mem = CVodeCreate(CV_BDF, CV_NEWTON);

    //CVodeInit allocates and initializes memory for a problem
    flag = CVodeInit(cvodes_mem, CvodesAlgorithm::RHSFunction, 0.0 /*initial time*/, branch_cvodes->y_);
    assert(flag==CV_SUCCESS);

    //specify integration tolerances. MUST be called before CVode.
    flag = CVodeSVtolerances(cvodes_mem, kRelativeTolerance, absolute_tolerance);
    assert(flag==CV_SUCCESS);

    //specify this branch as user data parameter to be past to functions
    flag = CVodeSetUserData(cvodes_mem, local);
    assert(flag==CV_SUCCESS);

    //specify g as root function and roots
#ifdef NDEBUG
    int roots_count = 1; //AP threshold
    int roots_direction [roots_count] = {1}; //root [0] only increasing
#else
    int roots_count = 3; //AP threshold + alarm for impossible minimum+max voltage
    int roots_direction[roots_count] = {1,0,0};
    //root [0] only increasing, [1] and [2] both ways
#endif
    flag = CVodeRootInit(cvodes_mem, roots_count, CvodesAlgorithm::RootFunction);
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
    void * cvodes_mem = branch_cvodes->cvodes_mem_;

#ifdef NDEBUG
    int roots_count = 1; //AP threshold
#else
    int roots_count = 3; //AP threshold + alarm for impossible minimum+max voltage
#endif
    int roots_found[roots_count];

    int flag=0;
    realtype tout = 0;
    hpx_t spikes_lco = HPX_NULL;

    while(local->nt_->_t < input_params_->tstop_)
    {
      //Beginning of execution of new event, update status
      CvodesAlgorithm::UpdateNrnThreadFromCvodeState(local);
      CvodesAlgorithm::ReevaluateBranch(local);

      //get tout as time of next undelivered event (if any)
      hpx_lco_sema_p(local->events_queue_mutex_);
      if (!local->events_queue_.empty())
          tout = local->events_queue_.top().first;
      hpx_lco_sema_v_sync(local->events_queue_mutex_);
      tout = std::min(input_params_->tstop_, tout);

      //call CVODE method: steps until reaching/passing tout;
      //TODO can it walk backwards for missed event?
      flag = CVode(branch_cvodes->cvodes_mem_, tout,
                   branch_cvodes->y_, &local->nt_->_t, CV_NORMAL);

      printf("At t = %0.4e   V =%14.6e  %14.6e  %14.6e\n",
             local->nt_->_t, NV_Ith_S(branch_cvodes->y_,0),
             NV_Ith_S(branch_cvodes->y_,1), NV_Ith_S(branch_cvodes->y_,2));

      if(flag==CV_ROOT_RETURN) //CVODE succeeded and roots found
      {
       flag = CVodeGetRootInfo(cvodes_mem, roots_found);
       assert(flag==CV_SUCCESS);
       assert(roots_found[0]!=0); //only possible root: AP threshold reached
#ifndef NDEBUG
       assert(roots_found[1]==0); //root can't be found or V too high
       assert(roots_found[2]==0); //root can't be found or V too low
#endif
       //NOTE: if root is found, integrator time is now at time of root!
       //(+1 value ascending, -1 valued descending)
       if (roots_found[0] > 0) //AP threshold reached from below
       {
           //TODO we can use a root to stop progress in time and wait for deliver?
           //use several hpx-all reduce to mark performance?
           spikes_lco=local->soma_->SendSpikes(local->nt_->_t);
       }
      }

      //Success from CVodeGetRootInfo or CVode
      if (flag == CV_SUCCESS)
        local->nt_->_t = tout;

      long num_steps=-1, num_jacob_evals=-1, num_rhs_evals=-1;
#ifdef NEUROX_CVODES_JACOBIAN_SOLVER==0
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
