#pragma once
#include "neurox.h"

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

namespace interpolators {

  class VariableTimeStep {
   public:
    VariableTimeStep();
    ~VariableTimeStep();

    /// absolute tolerance per equation
    N_Vector absolute_tolerance_;

    /// Initial values (voltages)
    N_Vector y_;

    /// CVODES structure
    void *cvodes_mem_;

    /// number of equations/vars in the system of ODEs
    /// i.e. compartments * equations per compartment
    int equations_count_;

    /// HPX actions registration
    static void RegisterHpxActions();

    /// mapping of y in CVODES to NrnThread->data
    double **state_var_map_;

    /// mapping of y in CVODES to NrnThread->data
    double **state_dv_map_;

    ///> Information of no-capacitance nodes
    class NoCapacitor
    {
      public:
        NoCapacitor() = delete;
        NoCapacitor(const Branch*);
        ~NoCapacitor();

        int * node_ids_; ///> no-cap node ids
        int * child_ids_; ///> id of nodes with no-cap parents
        int node_count_; ///> size of node_ids_
        int child_count_; ///> size of child_ids_
        Memb_list * memb_list_; ///> Memb_list of no-cap nodes
    } * no_cap_; ///> info on non-capacitors nodes

    static hpx_action_t Init;
    static hpx_action_t Run;
    static hpx_action_t Clear;

   private:

  /// CVODES BDF max-order
  const static int kBDFMaxOrder = 5;

  /// CVODES Mininum step size allowed (dt=0.025)
  constexpr static double kMinStepSize = 0; //13-6;

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
  static void GatherY(Branch *branch, N_Vector y);

  /// update ydot CVODES from NrnThread->data
  static void ScatterYdot(Branch *branch, N_Vector ydot);
  static void GatherYdot(Branch *branch, N_Vector ydot);

  static int RHSFunction(floble_t t, N_Vector y_, N_Vector ydot,
                         void *user_data);

  static void NoCapacitanceV(Branch * branch);

  /// g root function to compute g_i(t,y)
  static int RootFunction(realtype t, N_Vector y_, realtype *gout,
                          void *user_data);

  /// jacobian: compute J(t,y) on a dense matrix
  static int JacobianDense(long int N, floble_t t, N_Vector y_, N_Vector fy,
                              DlsMat J, void *user_data, N_Vector tmp1,
                              N_Vector tmp2, N_Vector tmp3);

  static int Init_handler();
  static int Run_handler();
  static int Clear_handler();
};

};  // interpolators

};  // neurox
