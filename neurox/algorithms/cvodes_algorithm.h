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

    /// mapping of y in CVODES to NrnThread->data
    double **state_var_map_;

    /// mapping of y in CVODES to NrnThread->data
    double **state_dv_map_;

    /// hpx for last sent spikes
    hpx_t spikes_lco_;

    /// HPX actions registration
    static void RegisterHpxActions();

    static hpx_action_t Init;
    static hpx_action_t Run;
    static hpx_action_t Clear;

    /// Data structure past as user-data argument to CVODES functions
    class UserData {
     public:
      UserData() = delete;
      UserData(Branch *);
      ~UserData();

      /// Branch this user data belongs to;
      Branch *branch_;

      /// VEC_D of last RHS call, before hines solver
      double *jacob_d_;

      /// temporary placeholder for data
      double *data_bak_;

      /// execution time of last RHS function call
      realtype rhs_last_time_;
      realtype rhs_second_last_time_;

    } * user_data_;

   private:
    static int Init_handler();
    static int Run_handler();
    static int Clear_handler();
  };

 private:
  /// CVODES Mininum step size allowed
  constexpr static double kMinStepSize = 1e-12;

  /// CVODES Relative torelance
  constexpr static double kRelativeTolerance = 1e-3;

  /// CVODES Absolute tolerance for voltage values
  constexpr static double kAbsToleranceVoltage = 1e-3;

  /// CVODES Absolute tolerance for mechanism states values
  constexpr static double kAbsToleranceMechStates = 1e-2;

  /// Time-window size for grouping of events to be delivered
  /// simmultaneously (0 for no grouping)
  constexpr static double kEventsDeliveryTimeWindow = 0.125;

  /// update NrnThread->data from with new CVODES state
  static void CopyYToVoltage(N_Vector y, Branch *branch);

  /// update CVODES from NrnThread->data
  static void CopyRHSToYdot(Branch *branch, N_Vector ydot);

  /// function defining the right-hand side function in y' = f(t,y).
  static int RHSFunction(floble_t t, N_Vector y_, N_Vector ydot,
                         void *user_data);

  /// g root function to compute g_i(t,y)
  static int RootFunction(realtype t, N_Vector y_, realtype *gout,
                          void *user_data);

  /// jacobian: compute J(t,y)
  static int JacobianFunction(long int N, floble_t t, N_Vector y_, N_Vector fy,
                              DlsMat J, void *user_data, N_Vector tmp1,
                              N_Vector tmp2, N_Vector tmp3);
};

};  // algorithm

};  // neurox
