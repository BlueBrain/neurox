#include "neurox/algorithms/cvodes_algorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

CvodesAlgorithm::CvodesAlgorithm():
data_(nullptr), absolute_tolerance_(nullptr)
{
}

CvodesAlgorithm::~CvodesAlgorithm() {}

const AlgorithmType CvodesAlgorithm::GetType() {
  return AlgorithmType::kCvodes;
}

const char* CvodesAlgorithm::GetTypeString() {
  return "CVODES";
}


/// f routine to compute y' = f(t,y).
int CvodesAlgorithm::F(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    realtype y1, y2, y3, yd1, yd3;

    y1 = NV_Ith_S(y,0); y2 = NV_Ith_S(y,1); y3 = NV_Ith_S(y,2);

    yd1 = NV_Ith_S(ydot,0) = RCONST(-0.04)*y1 + RCONST(1.0e4)*y2*y3;
    yd3 = NV_Ith_S(ydot,2) = RCONST(3.0e7)*y2*y2;
          NV_Ith_S(ydot,1) = -yd1 - yd3;

    return(0);
}

/// g routing to compute g_i(t,y) for i = 0,1.
int CvodesAlgorithm::G(realtype t, N_Vector y, realtype *gout, void *user_data)
{
    realtype y1, y3;

    y1 = NV_Ith_S(y,0); y3 = NV_Ith_S(y,2);
    gout[0] = y1 - RCONST(0.0001);
    gout[1] = y3 - RCONST(0.01);

    return(0);
}

//jacobian routine: compute J(t,y) = df/dy
int CvodesAlgorithm::Jacobian(long int N, realtype t,
               N_Vector y, N_Vector fy,
               DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    realtype y1, y2, y3;

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

/// resets the integrator at time t (NEURON book page 173)
void CvodesAlgorithm::Initialize(Branch * b, floble_t t)
{
    return;
}

/// performs an integration step to a new time t, and returns t (NEURON book page 173)
floble_t CvodesAlgorithm::Advance(Branch *b)
{
    return 0;
}

/// integration step that returns before next discontinuity event (NEURON book page 173)
floble_t CvodesAlgorithm::Interpolate(Branch *b)
{
    return 0;
}

void CvodesAlgorithm::Init() {

    int flag = CV_SUCCESS;

    //set-up data
    data_ = (UserData) malloc(sizeof *data_);
    data_->p[0]=3;

    y_ = N_VNew_Serial(kNumEquations);
    for (int i=0; i<kNumEquations; i++)
         NV_Ith_S(y_,i) = 0;

    //absolute tolerance array
    absolute_tolerance_ = N_VNew_Serial(kNumEquations);
    for (int i=0; i<kNumEquations; i++) //TODO set values
        NV_Ith_S(absolute_tolerance_, i) = kRelativeTolerance;

    //CVodeCreate creates an internal memory block for a problem to
    //be solved by CVODES, with Backward Differentiation (or Adams)
    //and Newton solver (recommended for stiff problems, see header)
    cvode_mem_ = CVodeCreate(CV_BDF, CV_NEWTON);

    //CVodeInit allocates and initializes memory for a problem
    flag = CVodeInit(cvode_mem_, CvodesAlgorithm::F, 0.0 /*initial time*/, y_);
    assert(flag==CV_SUCCESS);

    //specify integration tolerances. MUST be called before CVode.
    flag = CVodeSVtolerances(cvode_mem_, kRelativeTolerance, absolute_tolerance_);
    assert(flag==CV_SUCCESS);

    //specify g as root function with 2 components
    flag = CVodeRootInit(cvode_mem_, 2, CvodesAlgorithm::G);
    assert(flag==CV_SUCCESS);

    //initializes the memory record and sets various function
    //fields specific to the dense linear solver module.
    //TODO this is newton solver, shall we use direct solve instead?
    //Note: direct solvers give the solution (LU-decomposition, etc)
    //Indirect solvers require iterations (eg Jacobi method)
    flag = CVDense(cvode_mem_, kNumEquations);
    assert(flag==CV_SUCCESS);

    //specify the dense (user-supplied) Jacobian function. Compute J(t,y).
    flag = CVDlsSetDenseJacFn(cvode_mem_, CvodesAlgorithm::Jacobian);
    assert(flag==CV_SUCCESS);
}

void CvodesAlgorithm::Clear() {
    /* Free y vector */
    N_VDestroy_Serial(y_);

    /* Free integrator memory */
    CVodeFree(&cvode_mem_);
}

double CvodesAlgorithm::Launch() {
    realtype first_output_t = 0.4;
    int rootsfound[2];
    int flag=0;
    int iteration=0;
    while(1) {
      // call CVODE method
      // the CV_NORMAL task is to have the solver take internal steps until
      // it has reached or just passed the user specified tout parameter.
      flag = CVode(cvode_mem_, first_output_t, y_, &t_, CV_NORMAL);
      printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n",
             t, NV_Ith_S(y_,0), NV_Ith_S(y_,1), NV_Ith_S(y_,2));

      if(flag==CV_ROOT_RETURN)
      {
        //CVode succeeded, and found one or more roots.
        //If nrtfn > 1, call CVodeGetRootInfo to see
        //which g_i were found to have a root at (*tret).

       flag = CVodeGetRootInfo(cvode_mem_, rootsfound);
       assert(flag==CV_SUCCESS);
       printf(" rootsfound[0]=%d; rootsfound[1]=%d;\n", rootsfound[0], rootsfound[1]);

      }

      if (flag>0) break;

      if (flag == CV_SUCCESS) {
        iteration++;
        first_output_t *= 10; //output time factor
      }

      if (iteration == 12) break; //number of output times
    }

    /* Print some final statistics */
    //PrintFinalStats(cvode_mem_);
}

void CvodesAlgorithm::StepBegin(Branch*) {}

void CvodesAlgorithm::StepEnd(Branch* b, hpx_t spikesLco) {
}

void CvodesAlgorithm::Run(Branch* b, const void* args)
{
    
}

hpx_t CvodesAlgorithm::SendSpikes(Neuron* b, double tt, double) {
}
