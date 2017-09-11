#include "neurox/algorithms/cvodes_algorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

CvodesAlgorithm::CvodesAlgorithm():
data(nullptr)
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
int CvodesAlgorithm::f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype y1, y2, y3, yd1, yd3;
  UserData data;
  realtype p1, p2, p3;

  y1 = NV_Ith_S(y,0); y2 = NV_Ith_S(y,1); y3 = NV_Ith_S(y,2);
  data = (UserData) user_data;
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  yd1 = NV_Ith_S(ydot,0) = -p1*y1 + p2*y2*y3;
  yd3 = NV_Ith_S(ydot,2) = p3*y2*y2;
        NV_Ith_S(ydot,1) = -yd1 - yd3;

  return 0;
}

//jacobian routine: compute J(t,y).
int CvodesAlgorithm::jacobian(long int N, realtype t,
               N_Vector y, N_Vector fy,
               DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y1, y2, y3;
  UserData data;
  realtype p1, p2, p3;

  y1 = NV_Ith_S(y,0);
  y2 = NV_Ith_S(y,1);
  y3 = NV_Ith_S(y,2);
  data = (UserData) user_data;
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  DENSE_ELEM(J,0,0) = -p1;
  DENSE_ELEM(J,0,1) = p2*y3;
  DENSE_ELEM(J,0,2) = p2*y2;
  DENSE_ELEM(J,1,0) =  p1;
  DENSE_ELEM(J,1,1) = -p2*y3-2*p3*y2;
  DENSE_ELEM(J,1,2) = -p2*y2;
  //todo wheres 2,0 and 2,2?
  DENSE_ELEM(J,2,1) = 2*p3*y2;

  return(0);
}

// fQ is the user-provided integrand routine. Compute fQ(t,y).
int CvodesAlgorithm::fQ(realtype t, N_Vector y, N_Vector qdot, void *user_data)
{
  NV_Ith_S(qdot,0) = NV_Ith_S(y,2);
  return 0;
}

void CvodesAlgorithm::Init() {

    int flag = CV_SUCCESS;

    //set-up data
    data = (UserData) malloc(sizeof *data);
    data->p[0]=3;

    //y0 is the initial condition vector y(t0)
    y0 = N_VNew_Serial(kNumEquations);
    for (int i=0; i<kNumEquations; i++)
         NV_Ith_S(y0,i) = 0;

    //yQ0 is initial condition of the quadrature values
    yQ0 = N_VNew_Serial(1);
    NV_Ith_S(yQ0,0) = 0.0;

    //t0 is initial time
    floble_t t0 = 0.0;

    //CVodeCreate creates an internal memory block for a problem to
    //be solved by CVODES, with Backward Differentiation (or Adams)
    //and Newton solver (recommended for stiff problems, see header)
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

    //CVodeInit allocates and initializes memory for a problem
    flag = CVodeInit(cvode_mem, CvodesAlgorithm::f, t0, y0);
    assert(flag==CV_SUCCESS);

    //specify integration tolerances. MUST be called before CVode.
    //TODO use CVodeWFtolerances to set diff absolute tolerances per term?
    flag = CVodeSStolerances(cvode_mem, kRelativeTolerance, kAbsoluteTolerance);
    assert(flag==CV_SUCCESS);

    //set a pointer to user data that will be passed to the user's f
    //function every time f is called.
    flag = CVodeSetUserData(cvode_mem, data);
    assert(flag==CV_SUCCESS);

    //initializes the memory record and sets various function
    //fields specific to the dense linear solver module.
    //TODO we should provide our own solver here, see page 10 of CVODES guide!
    flag = CVDense(cvode_mem, kNumEquations);
    assert(flag==CV_SUCCESS);

    //specify the dense Jacobian function. Compute J(t,y).
    flag = CVDlsSetDenseJacFn(cvode_mem, CvodesAlgorithm::jacobian);
    assert(flag==CV_SUCCESS);

    //allocates and initializes quadrature related memory for a problem
    //(quadrature = area = integral in finit or infinite interval)
    flag = CVodeQuadInit(cvode_mem, fQ, yQ0);
    assert(flag==CV_SUCCESS);

    //sets quadrature tolerance
    flag = CVodeQuadSStolerances(cvode_mem, kRelativeTolerance, kAbsoluteToleranceQuadrature);
    assert(flag==CV_SUCCESS);

    //Quadrature optional input specification functions:
    //are quadrature variables considered in the error control?
    flag = CVodeSetQuadErrCon(cvode_mem, TRUE);
    assert(flag==CV_SUCCESS);

}

void CvodesAlgorithm::Clear() {
}

double CvodesAlgorithm::Launch() {
}

void CvodesAlgorithm::StepBegin(Branch*) {}

void CvodesAlgorithm::StepEnd(Branch* b, hpx_t spikesLco) {
}

void CvodesAlgorithm::Run(Branch* b, const void* args)
{
    
}

hpx_t CvodesAlgorithm::SendSpikes(Neuron* b, double tt, double) {
}
