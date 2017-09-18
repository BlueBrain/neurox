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

  const AlgorithmId GetId() override;
  const char* GetString() override;

  void Init() override;
  void Clear() override;
  double Launch() override;

  void StepBegin(Branch*) override;
  void StepEnd(Branch*, hpx_t) override;
  void Run(Branch*, const void*) override;
  hpx_t SendSpikes(Neuron*, double, double) override;

  class BranchCvodes : public AlgorithmMetadata
  {
    public:

      BranchCvodes();
      ~BranchCvodes();

      /// absolute tolerance per equation
      N_Vector absolute_tolerance_;

      /// Initial values (voltages)
      N_Vector y_;

      /// CVODES structure
      void *cvodes_mem_;

      /// number of equations/vars in the system of ODEs
      /// i.e. compartments * equations per compartment
      int equations_count_;

      /// number of states per compartment
      int equations_per_compartment_;

      /// mapping of equations y to NrnThread->data
      /// (similar to ode_map in NEURON)
      double **equations_map_;

      /// minimum step size
      static double min_step_size_;

      /// HPX actions registration
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
  constexpr static double kRelativeTolerance = 1e-3;
  constexpr static double kAbsToleranceVoltage = 1e-3;
  constexpr static double kAbsToleranceMechStates = 1e-2;

  /// function that reavaluates all elementes in NrnThread->data
  static int ReevaluateBranch(Branch * branch);

  /// function defining the right-hand side function in y' = f(t,y).
  static int RHSFunction(floble_t t, N_Vector y_, N_Vector ydot, void *user_data);

  /// g root function to compute g_i(t,y)
  static int RootFunction(realtype t, N_Vector y_, realtype *gout, void *user_data);

  /// jacobian: compute J(t,y)
  static int JacobianSparseMatrix(
          realtype t,
          N_Vector y, N_Vector fy, SlsMat JacMat, void *user_data,
          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

  /// jacobian: compute J(t,y)
  static int JacobianFunction(long int N, floble_t t,
                      N_Vector y_, N_Vector fy,
                      DlsMat J, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

};

};  // algorithm

};  // neurox
