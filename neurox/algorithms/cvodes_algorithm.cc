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


/// f routine to compute y' = f(t,y).
int CvodesAlgorithm::RHSFunction(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    Branch * branch = (Branch*) user_data;

    branch->CallModFunction(Mechanism::ModFunctions::kCurrent);
    //TODO this is an exception
    //b->CallModFunction(Mechanism::ModFunctions::kCurrentCapacitance);

   // double a = _vec_shadow[i]; //computed by current functions

    realtype y1, y2, y3, yd1, yd3;

    y1 = NV_Ith_S(y,0); y2 = NV_Ith_S(y,1); y3 = NV_Ith_S(y,2);

    yd1 = NV_Ith_S(ydot,0) = RCONST(-0.04)*y1 + RCONST(1.0e4)*y2*y3;
    yd3 = NV_Ith_S(ydot,2) = RCONST(3.0e7)*y2*y2;
          NV_Ith_S(ydot,1) = -yd1 - yd3;

    return(0);
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

double get_g()
{
    return 0.0;
}

//jacobian routine: compute J(t,y) = df/dy
int CvodesAlgorithm::JacobianDenseFunction(
        long int N, realtype t,
        N_Vector y, N_Vector fy,
        DlsMat J, void *user_data,
        N_Vector, N_Vector, N_Vector)
{
    realtype y1, y2, y3;
    Branch * branch = (Branch*) user_data;
    realtype * ydata = N_VGetArrayPointer_Serial(y);
    assert(ydata==branch->nt_->_data);
    assert(N==branch->nt_->_ndata);
    assert(t==branch->nt_->_t);

    //Jacobian for nt->data includes
    //dV_0/dt, ..., dV_n/dt, i_0, i_1, ..., i_n
    //constants have jacobi=0
    //currents have jacobi=g (given by nrn_jacob)

    //CVODES guide: must load NxN matrix J with the approximation
    //of Jacobian J(t,y) at point (t,y)


    realtype ** jacob = J->cols;
    int v_offset = branch->nt_->_actual_v - branch->nt_->_data;
    realtype ** jacob_v = &jacob[v_offset];

    //TODO we can just set 'dt' to 1
    //Scale derivatives to a one time-unit size
    const realtype rev_dt = input_params_->rev_dt_;

    //add parent/children compartments currents to C*dV/dt
    //add rhs: includes all other currents
    int compartments_count = branch->nt_->end;
    realtype *a = branch->nt_->_actual_a ;
    realtype *b = branch->nt_->_actual_b ;
    realtype *rhs = branch->nt_->_actual_rhs ;

    for (int c=0; c<compartments_count; c++)
       jacob_v[c][0] += a[c]*rev_dt;

    for (int c=0; c<compartments_count; c++)
       jacob_v[c][0] += b[c]*rev_dt;

    for (int c=0; c<compartments_count; c++)
       jacob_v[c][0] += rhs[c]*rev_dt;

    int compartment_id=-1;
    int g_index=-1;
    Memb_list * mech_instances= nullptr;
    for (int m=0; m<mechanisms_count_; m++)
    {
        mech_instances = &branch->mechs_instances_[m];
        for (int i=0; i<mech_instances->nodecount; b++)
        {
          compartment_id=mech_instances->nodeindices[i];

          //updates the main current functions dV/dt
          jacob_v[compartment_id][i] += get_g()*rev_dt; // g == di/dV

          //update di/dV

          //update rhs?
        }
    }

    y1 = NV_Ith_S(y,0); y2 = NV_Ith_S(y,1); y3 = NV_Ith_S(y,3);
    DENSE_ELEM(J,0,0) = RCONST(-0.04);
    DENSE_ELEM(J,0,1) = RCONST(1.0e4)*y3;
    DENSE_ELEM(J,0,2) = RCONST(1.0e4)*y2;
    DENSE_ELEM(J,1,0) = RCONST(0.04);
    DENSE_ELEM(J,1,1) = RCONST(-1.0e4)*y3-RCONST(6.0e7)*y2;
    DENSE_ELEM(J,1,2) = RCONST(-1.0e4)*y2;
    DENSE_ELEM(J,2,1) = RCONST(6.0e7)*y2;

    return(0);
}

/////////////////////// Algorithm abstract class ////////////////////////


void CvodesAlgorithm::Init() {
  BranchCvodes::min_step_size_ = input_params_->dt_;
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


void CvodesAlgorithm::StepBegin(Branch*)
{

}

void CvodesAlgorithm::StepEnd(Branch* b, hpx_t spikesLco) {
}

void CvodesAlgorithm::Run(Branch* b, const void* args)
{

}

hpx_t CvodesAlgorithm::SendSpikes(Neuron* n, double tt, double) {
    return Neuron::SendSpikesAsync(n, tt);
}


//////////////////////////// BranchCvodes /////////////////////////


CvodesAlgorithm::BranchCvodes::BranchCvodes()
    :cvodes_mem_(nullptr), iterations_count_(0)
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
    N_Vector y = branch_cvodes->y_;
    void * cvode_mem = branch_cvodes->cvodes_mem_;

    int flag = CV_ERR_FAILURE;

    //equations: one per vdata (some will be constants, with jacobian=0)
    int num_equations = local->nt_->_ndata;

    //create serial vector for y
    //TODO guide page 157, we can have mem-protected N_Vectors
    y = N_VMake_Serial(num_equations, local->nt_->_data);

    //absolute tolerance array
    absolute_tolerance = N_VNew_Serial(num_equations);
    for (int i=0; i<num_equations; i++) //TODO set values
        NV_Ith_S(absolute_tolerance, i) = kRelativeTolerance;

    //CVodeCreate creates an internal memory block for a problem to
    //be solved by CVODES, with Backward Differentiation (or Adams)
    //and Newton solver (recommended for stiff problems, see header)
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

    //CVodeInit allocates and initializes memory for a problem
    flag = CVodeInit(cvode_mem, CvodesAlgorithm::RHSFunction, 0.0 /*initial time*/, y);
    assert(flag==CV_SUCCESS);

    //specify integration tolerances. MUST be called before CVode.
    flag = CVodeSVtolerances(cvode_mem, kRelativeTolerance, absolute_tolerance);
    assert(flag==CV_SUCCESS);

    //specify this branch as user data parameter to be past to functions
    flag = CVodeSetUserData(cvode_mem, local);
    assert(flag==CV_SUCCESS);

    //specify g as root function and roots
#ifdef NDEBUG
    int roots_count = 1; //AP threshold
#else
    int roots_count = 3; //AP threshold + alarm for impossible minimum+max voltage
#endif
    flag = CVodeRootInit(cvode_mem, roots_count, CvodesAlgorithm::RootFunction);
    assert(flag==CV_SUCCESS);

    //initializes the memory record and sets various function
    //fields specific to the dense linear solver module.
    //TODO this is newton solver, shall we use direct solve instead?
    //Note: direct solvers give the solution (LU-decomposition, etc)
    //Indirect solvers require iterations (eg Jacobi method)
    //sundials_direct.h wraps solvers??

    if (kSparseMatrix)
    {
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
        assert(0);
    }
    else
    {
      flag = CVDense(cvode_mem, num_equations);
      assert(flag==CV_SUCCESS);

      //specify the dense (user-supplied) Jacobian function. Compute J(t,y).
      flag = CVDlsSetDenseJacFn(cvode_mem, CvodesAlgorithm::JacobianDenseFunction);
      assert(flag==CV_SUCCESS);
    }

    //TODO added Bruno (see cvs_guide.pdf page 43)
    CVodeSetInitStep(cvode_mem, min_step_size_);
    CVodeSetMinStep(cvode_mem, min_step_size_);
    CVodeSetMaxStep(cvode_mem, CoreneuronAlgorithm::CommunicationBarrier::kCommStepSize);
    CVodeSetStopTime(cvode_mem, input_params_->tstop_);
    CVodeSetMaxOrd(cvode_mem, 5); //max order of the BDF method

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
    TimedEvent * top_event = nullptr;

    while(local->nt_->_t < input_params_->tstop_)
    {
      //get tout i.e. max time of next step
      hpx_lco_sema_p(local->events_queue_mutex_);
      if (!local->events_queue_.empty())
          tout = local->events_queue_.top().first;
      hpx_lco_sema_v_sync(local->events_queue_mutex_);

      //time of next step
      tout = std::min(input_params_->tstop_, tout);

      //call CVODE method: take internal steps until it has
      //reached or just passed tout; update t;
      //TODO can it walk backwards for missed event?
      flag = CVode(branch_cvodes->cvodes_mem_, tout,
                   branch_cvodes->y_, &local->nt_->_t, CV_NORMAL);

      printf("At t = %0.4e   V =%14.6e  %14.6e  %14.6e\n",
             local->nt_->_t, NV_Ith_S(branch_cvodes->y_,0),
             NV_Ith_S(branch_cvodes->y_,1), NV_Ith_S(branch_cvodes->y_,2));

      if(flag==CV_ROOT_RETURN) //CVODE succeeded and roots found
      {
        //CVode succeeded, and found roots
        //(+1 value ascending, -1 valued descending)
       flag = CVodeGetRootInfo(cvodes_mem, roots_found);
       assert(flag==CV_SUCCESS);
       assert(roots_found[0]!=0); //only possible root: AP threshold reached
#ifndef NDEBUG
       assert(roots_found[1]==0); //root can't be found or V too high
       assert(roots_found[2]==0); //root can't be found or V too low
#endif
       //if root is found, integrator time is at time of root!
       if (roots_found[0] > 0) //AP threshold reached from below
       {
           //TODO
           //use several hpx-all reduce to mark performance?
           spikes_lco=local->soma_->SendSpikes(local->nt_->_t);
       }
      }
      else if (flag == CV_SUCCESS)
        branch_cvodes->iterations_count_++;

      //delivers all events whithin the next min step size
      local->DeliverEvents(local->nt_->_t + min_step_size_);

      //No discontinuities: only sets vecplay->*pd for next discont. event
      local->FixedPlayContinuous();
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

    N_VDestroy_Serial(branch_cvodes->y_);  /* Free y vector */
    CVodeFree(&branch_cvodes->cvodes_mem_); /* Free integrator memory */
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
