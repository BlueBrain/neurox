#pragma once
#include "neurox.h"

//by cvodes/examples/cvsRoberts_ASAi_dns.c
#include "cvodes/cvodes.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_types.h"
#include "sundials/sundials_math.h"
#include "cvodes/cvodes_dense.h"

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

 private:
  static const int kNumEquations = 3;
  constexpr static double kRelativeTolerance = 1e-3;
  constexpr static double kAbsoluteTolerance = 1e-3;
  constexpr static double kAbsoluteToleranceQuadrature = 1e-3;

  typedef struct {
    floble_t p[kNumEquations];
  } *UserData;

  UserData data;

  /// Initial condition
  N_Vector y0;

  /// initial condition of the quadrature values
  N_Vector yQ0;

  void *cvode_mem;

  /// function defining the right-hand side function in y' = f(t,y).
  static int f(floble_t t, N_Vector y0, N_Vector ydot, void *user_data);

  /// jacobian: compute J(t,y)
  static int jacobian(long int N, floble_t t,
                      N_Vector y, N_Vector fy,
                      DlsMat J, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

  /// fQ is the user-provided integrand routine. Compute fQ(t,y).
  static int fQ(realtype t, N_Vector y, N_Vector qdot, void *user_data);
};

};  // algorithm

};  // neurox
