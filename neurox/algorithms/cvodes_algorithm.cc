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
int CvodesAlgorithm::RHSFunction(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    Branch * branch = (Branch*) user_data;
    BranchCvodes* branch_cvodes = (BranchCvodes*) branch->soma_->algorithm_metadata_;
    assert(NV_LENGTH_S(y) == NV_LENGTH_S(ydot));


    //changes have to be made in ydot, y is the state at the previous step
    //so we copy y to ydot, and apply the changes to ydot
    memcpy(NV_DATA_S(ydot), NV_DATA_S(y), NV_LENGTH_S(y)*sizeof(floble_t));

    //update vecplay pointers to point to right place
    VecplayContinuousX * vecplay = nullptr;
    int vecplay_pd_offset=-1;
    for (int v=0; v<branch->nt_->n_vecplay; v++)
    {
        vecplay=(VecplayContinuousX*) branch->nt_->_vecplay[v];
        vecplay_pd_offset = vecplay->pd_ - branch->nt_->_data;
        vecplay->pd_ = &NV_DATA_S(ydot)[vecplay_pd_offset];
    }
    branch->nt_->_data = NV_DATA_S(ydot);

    //set time step to the time we want to jump to
    assert(t > branch->nt_->_t);
    branch->nt_->_dt = t - branch->nt_->_t;
    NV_DATA_S(ydot)[branch_cvodes->equations_count_-1] = t;
    assert(branch->nt_->_dt >= BranchCvodes::min_step_size_);

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

    return 0;
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

//jacobian routine: compute J(t,y) = df/dy
int CvodesAlgorithm::JacobianFunction(
        long int N, realtype t,
        N_Vector y, N_Vector fy,
        DlsMat J, void *user_data,
        N_Vector, N_Vector, N_Vector)
{

    Branch * branch = (Branch*) user_data;
    BranchCvodes* branch_cvodes = (BranchCvodes*) branch->soma_->algorithm_metadata_;

    int a_offset = branch->nt_->_actual_a - branch->nt_->_data;
    int b_offset = branch->nt_->_actual_b - branch->nt_->_data;
    int v_offset = branch->nt_->_actual_v - branch->nt_->_data;
    int d_offset = branch->nt_->_actual_d - branch->nt_->_data;
    int rhs_offset = branch->nt_->_actual_rhs - branch->nt_->_data;
    int compartments_count = branch->nt_->end;

    //TODO double-check: Jacob has been stored in f(y,y), not y???
    realtype *y_data = N_VGetArrayPointer_Serial(y);
    realtype *fy_data = N_VGetArrayPointer_Serial(fy);
    realtype *a = &y_data[a_offset];
    realtype *b = &y_data[b_offset];
    realtype *v = &y_data[v_offset];
    realtype *d = &y_data[d_offset];
    realtype *rhs = &y_data[rhs_offset];
    int *p = branch->nt_->_v_parent_index;
    assert(y_data==branch->nt_->_data);
    assert(t==branch->nt_->_t);

    //Jacobian for nt->data includes all voltages, currents,
    //ions states and point-processes wreights

    //CVODES guide: must load NxN matrix J with the approximation
    //of Jacobian J(t,y) at point (t,y)

    realtype ** jacob = J->cols;
    realtype ** jacob_d = &jacob[d_offset];
    realtype ** jacob_rhs = &jacob[rhs_offset];

    //Scale factor for derivatives (based on previous step taken)
    const int t_index = branch_cvodes->equations_count_-1;
    assert(fy_data[t_index] > y_data[t_index]);
    realtype dt = fy_data[t_index] - y_data[t_index];
    const realtype rev_dt = 1 / dt;


    //add constributions from parent/children compartments
    //TODO fix positions
    floble_t dv=-1;
    for (offset_t i = 1; i < compartments_count; i++)
    {
      //from HinesSolver::SetupMatrixRHS
      dv = v[p[i]] - v[i];
      jacob_rhs[i][i] -= b[i] * dv * rev_dt;
      jacob_rhs[p[i]][i] += a[i] * dv * rev_dt;

      //from HinesSolver::SetupMatrixDiagonal
      jacob_d[i][i] -= b[i] * rev_dt;
      jacob_d[p[i]][i] -= a[i] * rev_dt;
    }


    //Add contributions from mechanisms

    for (int c=0; c<compartments_count; c++)
       jacob_d[c][0] += a[c]*rev_dt;

    for (int c=0; c<compartments_count; c++)
       jacob_d[c][0] += b[c]*rev_dt;

    for (int c=0; c<compartments_count; c++)
       jacob_d[c][0] += rhs[c]*rev_dt;

    int compartment_id=-1, g_data_offset=-1;
    int data_offset=  tools::Vectorizer::SizeOf(compartments_count)*6;
    realtype g=-1;
    for (int m = 0; m < neurox::mechanisms_count_; m++)
    {
        Mechanism * mech=mechanisms_[m];
        Memb_list* mech_instances = &branch->mechs_instances_[m];

        //ions, insert updates from second order corrent
        if (mech->is_ion_)
        {
            //TODO not needed, its inserted below?

    }
        //if no current function, no current jacobian
        else if (mech->memb_func_.current != NULL)
        {
            //TODO add for capacitance
          //index of g is "typically" size of data -1;
          //TODO hard-coded expection
          int g_index = mech->data_size_-1;
          if (mech->type_ == MechanismTypes::kProbAMPANMDA_EMS
                || mech->type_ == MechanismTypes::kProbGABAAB_EMS
                || mech->type_ == MechanismTypes::kExpSyn)
              g_index = mech->data_size_-2;

          for (int n=0; n<mech_instances->nodecount; n++)
          {
#if LAYOUT == 1
            g_data_offset = data_offset + mech->data_size_ * n + g_index;
#else
            g_data_offset = data_offset + tools::Vectorizer::SizeOf(mech_instances->nodecount) * g_index + n;
#endif
            //updates the main current functions dV/dt
            g = y_data[data_offset];
            compartment_id = mech_instances->nodeindices[m];
            jacob_d[compartment_id][g_data_offset] += g*rev_dt; // g == di/dV
          }
        }
        data_offset = tools::Vectorizer::SizeOf(mech_instances->nodecount)*mech->data_size_;
    }


    //update of weights of VecPlayContinuous based on time
    VecplayContinuousX * vecplay = nullptr;
    int vecplay_pd_offset=-1;
    for (int v=0; v<branch->nt_->n_vecplay; v++)
    {
        vecplay=(VecplayContinuousX*) branch->nt_->_vecplay[v];
        vecplay_pd_offset = vecplay->pd_ - branch->nt_->_data;
        jacob[vecplay_pd_offset][t_index] =
                 ( vecplay->Interpolate(fy_data[t_index])
                  -vecplay->Interpolate(y_data[t_index])
                 ) / dt;
    }

    //jacobian for time
    jacob[t_index][t_index] = dt;

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
    :cvodes_mem_(nullptr), iterations_count_(-1), equations_count_(-1)
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
    int & equations_count = branch_cvodes->equations_count_;
    NrnThread *& nt = local->nt_;

    int flag = CV_ERR_FAILURE;

    //equations: one per data, weight and time (constant have jacob=0)
    equations_count = local->nt_->_ndata + local->nt_->n_weight + 1;
    branch_cvodes->iterations_count_=0;

    //create array y for state: nt->data, nt->weights and time
    floble_t * y_data = new floble_t[equations_count];
    std::copy(nt->_data, nt->_data+local->nt_->_ndata, y_data);
    std::copy(nt->weights, nt->weights + nt->n_weight, y_data + nt->_ndata);
    y_data[equations_count=-1]=local->nt_->_t;
    tools::Vectorizer::Delete(nt->_data);
    delete[] nt->weights;
    nt->_data = y_data;
    nt->weights = &y_data[nt->_ndata];
    branch_cvodes->y_ = N_VMake_Serial(equations_count, y_data);

    //absolute tolerance array
    absolute_tolerance = N_VNew_Serial(equations_count);
    for (int i=0; i<equations_count; i++) //TODO set values
        NV_Ith_S(absolute_tolerance, i) = kRelativeTolerance;

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
#else
    int roots_count = 3; //AP threshold + alarm for impossible minimum+max voltage
#endif
    flag = CVodeRootInit(cvodes_mem, roots_count, CvodesAlgorithm::RootFunction);
    assert(flag==CV_SUCCESS);

    //initializes the memory record and sets various function
    //fields specific to the dense linear solver module.
    //TODO this is newton solver, shall we use direct solve instead?
    //Note: direct solvers give the solution (LU-decomposition, etc)
    //Indirect solvers require iterations (eg Jacobi method)
    //sundials_direct.h wraps solvers??

    /* install superlumt and klu first

    // Call CVSuperLUMT to specify the CVSuperLUMT sparse direct linear solver
    int nnz = kNumEquations * kNumEquations;
    flag = CVSuperLUMT(cvode_mem, 1, kNumEquations, nnz);
    //TODO? flag = CVKLU(cvode_mem, 1, kNumEquations, nnz);
    assert(flag==CV_SUCCESS);

    //specify the dense (user-supplied) Jacobian function. Compute J(t,y).
    flag = CVSlsSetSparseJacFn(cvode_mem_, CvodesAlgorithm::JacobianSparseMatrix);
    assert(flag==CV_SUCCESS);
    */

    flag = CVDense(cvodes_mem, equations_count);
    assert(flag==CV_SUCCESS);

    //specify the dense (user-supplied) Jacobian function. Compute J(t,y).
    flag = CVDlsSetDenseJacFn(cvodes_mem, CvodesAlgorithm::JacobianFunction);
    assert(flag==CV_SUCCESS);

    CVodeSetInitStep(cvodes_mem, min_step_size_);
    CVodeSetMinStep(cvodes_mem, min_step_size_);
    CVodeSetMaxStep(cvodes_mem, CoreneuronAlgorithm::CommunicationBarrier::kCommStepSize);
    CVodeSetStopTime(cvodes_mem, input_params_->tstop_);
    CVodeSetMaxOrd(cvodes_mem, 5); //max order of the BDF method

    //see chapter 6.4 -- User supplied functions
    //see chapter 8 -- providing alternate linear solver modules

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
      //delivers all events whithin the next min step size
      local->DeliverEvents(local->nt_->_t + min_step_size_);

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
      {
        branch_cvodes->iterations_count_++;
        local->nt_->_t = tout;
      }
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
