#pragma once
#include "neurox.h"

// COPIED FROM cvsRoberts_dns.c (dense) and cvsRoberts_sps.c (sparse matrix)
// NO SENSITIVITY ANALYSIS: FSA or ASA avaiable at cvsRoberts_ASA_idns.c

#include <cvodes/cvodes.h>           /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h> /* definition of type realtype */

// For Sparse Matrix resolutions
//#include <cvodes/cvodes_superlumt.h>  /* prototype for CVSUPERLUMT */
#include <sundials/sundials_sparse.h> /* definitions SlsMat */

// For Dense Matrix resolutions
#include <cvodes/cvodes_dense.h>     /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */

// For Approx Diagonal matrix
#include <cvodes/cvodes_diag.h>

using namespace neurox;

namespace neurox {

namespace algorithms {

class CvodesAlgorithm : public Algorithm {
 public:
  CvodesAlgorithm();
  ~CvodesAlgorithm();

  const AlgorithmId GetId() override;
  const char *GetString() override;

  void Init() override;
  void Clear() override;
  double Launch() override;

  void StepBegin(Branch *) override;
  void StepEnd(Branch *, hpx_t) override;
  void Run(Branch *, const void *) override;
  hpx_t SendSpikes(Neuron *, double, double) override;

  class BranchCvodes : public AlgorithmMetadata {
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

    /// hpx for last sent spikes
    hpx_t spikes_lco_;

    /// HPX actions registration
    static void RegisterHpxActions();

    /// temporary placeholder for data
    double *data_bak_;

    /// mapping of y in CVODES to NrnThread->data
    double **state_var_map_;

    /// mapping of y in CVODES to NrnThread->data
    double **state_dv_map_;

    static hpx_action_t Init;
    static hpx_action_t Run;
    static hpx_action_t Clear;

   private:
    static int Init_handler();
    static int Run_handler();
    static int Clear_handler();
  };

 private:
  /// CVODES BDF max-order
  const static int kBDFMaxOrder = 5;

  /// CVODES Mininum step size allowed (dt=0.025)
  constexpr static double kMinStepSize = 1e-3;

  /// CVODES Relative torelance
  constexpr static double kRelativeTolerance = 1e-3;

  /// CVODES Absolute tolerance for voltage values
  constexpr static double kAbsToleranceVoltage = 1e-3;

  /// CVODES Absolute tolerance for mechanism states values
  constexpr static double kAbsToleranceMechStates = 1e-3;

  /// Time-window size for grouping of events to be delivered
  /// simmultaneously (0 for no grouping)
  constexpr static double kEventsDeliveryTimeWindow = 0.125;

  /// update NrnThread->data from with new y state
  static void ScatterY(Branch *branch, N_Vector y);

  /// update CVODES from NrnThread->data
  static void ScatterYdot(Branch *branch, N_Vector ydot);

  /// update CVODES from NrnThread->data
  static void GatherYdot(Branch *branch, N_Vector ydot);

  static int RHSFunction(floble_t t, N_Vector y_, N_Vector ydot,
                         void *user_data);

  /// g root function to compute g_i(t,y)
  static int RootFunction(realtype t, N_Vector y_, realtype *gout,
                          void *user_data);

  /// jacobian: compute J(t,y) on a dense matrix
  static int JacobianDense(long int N, floble_t t, N_Vector y_, N_Vector fy,
                              DlsMat J, void *user_data, N_Vector tmp1,
                              N_Vector tmp2, N_Vector tmp3);
};

};  // algorithm

};  // neurox
