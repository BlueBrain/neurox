#pragma once
#include "neurox.h"

//COPIED FROM cvsRoberts_dns.c (dense) and cvsRoberts_sps.c (sparse matrix)
//NO SENSITIVITY ANALYSIS: FSA or ASA avaiable at cvsRoberts_ASA_idns.c

#include <cvodes/cvodes.h>            /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>   /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h> /* definition of type realtype */

//For Sparse Matrix resolutions
//#include <cvodes/cvodes_superlumt.h>  /* prototype for CVSUPERLUMT */
#include <sundials/sundials_sparse.h> /* definitions SlsMat */

//For Dense Matrix resolutions
#include <cvodes/cvodes_dense.h>     /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */

using namespace neurox;

namespace neurox {

namespace algorithms {

class CvodesAlgorithm : public Algorithm {
 public:
  CvodesAlgorithm();
  ~CvodesAlgorithm();

  const AlgorithmType GetType() override;
  const char* GetTypeString() override;

  void Init() override;
  void Clear() override;
  double Launch() override;

  void StepBegin(Branch*) override;
  void StepEnd(Branch*, hpx_t) override;
  void Run(Branch*, const void*) override;
  hpx_t SendSpikes(Neuron*, double, double) override;

  class CvodesInfo : public AlgorithmMetadata {
   public:
    CvodesInfo();
    ~CvodesInfo();
  };

  class BranchCvodes
  {
    public:

      typedef struct {
        realtype *p; //size kNumEquations
      } *UserData;

      UserData data_;

      /// absolute tolerance per equation
      N_Vector absolute_tolerance_;

      /// Initial condition
      N_Vector y_;

      void *cvode_mem_;

      static void RegisterHpxActions();

      static hpx_action_t Init;
      static hpx_action_t Run;
      static hpx_action_t Clear;

    private:
      static int Init_handler();
      static int Run_handler();
      static int Clear_handler();
  };

 private:
  static const bool kSparseMatrix = false;
  static const int kNumEquations = 3;
  constexpr static double kRelativeTolerance = 1e-3;

  /// function defining the right-hand side function in y' = f(t,y).
  static int F(floble_t t, N_Vector y_, N_Vector ydot, void *user_data);

  /// g root function to compute g_i(t,y)
  static int G(realtype t, N_Vector y_, realtype *gout, void *user_data);

  /// jacobian: compute J(t,y)
  static int JacobianSparseMatrix(
          realtype t,
          N_Vector y, N_Vector fy, SlsMat JacMat, void *user_data,
          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

  /// jacobian: compute J(t,y)
  static int JacobianDenseMatrix(long int N, floble_t t,
                      N_Vector y_, N_Vector fy,
                      DlsMat J, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

  /// resets the integrator at time t (NEURON book page 173)
  static void Initialize(Branch * b, floble_t t);

  /// performs an integration step to a new time t, and returns t (NEURON book page 173)
  static floble_t Advance(Branch *b);

  /// integration step that returns before next discontinuity event (NEURON book page 173)
  static floble_t Interpolate(Branch *b);

};

};  // algorithm

};  // neurox
