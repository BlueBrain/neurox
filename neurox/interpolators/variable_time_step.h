#pragma once
#include "neurox.h"

#include <cvodes/cvodes.h>           /* prototypes for CVODE fcts, consts*/
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts, macros*/
#include <sundials/sundials_types.h> /* definition of type realtype*/
#include "cvodes/cvodes_impl.h"      /* definition of CVodeMem*/

// For Sparse Matrix resolutions
//#include <cvodes/cvodes_superlumt.h>  /* prototype for CVSUPERLUMT */
#include <sundials/sundials_sparse.h> /* definitions SlsMat */

// For Dense Matrix resolutions
#include <cvodes/cvodes_dense.h>     /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */

// For Pre-conditioned matrix solvers
#include <cvodes/cvodes_diag.h>      /*For Approx Diagonal matrix*/

using namespace neurox;

namespace neurox {

namespace interpolators {

class VariableTimeStep : public Interpolator {
 public:

  VariableTimeStep();
  ~VariableTimeStep();

  const char* GetString() override;
  void Init(Branch*)  override;
  void StepTo(Branch*, const double)  override;
  void Clear(Branch*) override;

  static void PrintStatistics(const Branch*);

  // Information of no-capacitance nodes
  int* no_cap_node_ids_;        ///> no-cap node ids
  int  no_cap_node_ids_count_;  ///> size of node_ids_
  int* no_cap_child_ids_;       ///> id of nodes with no-cap parents
  int  no_cap_child_ids_count_; ///> size of child_ids_
  Memb_list *no_cap_ml_;  ///> Memb_list of no-cap nodes

private:

  /// absolute tolerance per equation
  N_Vector absolute_tolerance_;

  /// Initial values (voltages)
  N_Vector y_;

  /// CVODES structure
  CVodeMem cvode_mem_;

  /// number of equations/vars in the system of ODEs
  /// i.e. compartments * equations per compartment
  int equations_count_;

  /// mapping of y in CVODES to NrnThread->data
  double **state_var_map_;

  /// mapping of y in CVODES to NrnThread->data
  double **state_dv_map_;

  /// possible operations between NrnThread->data and CVodes states
  typedef enum CopyOps {
    kScatterY,
    kGatherY,
    kScatterYdot,
    kGatherYdot
  } CopyOp;

  /// copy data to/from branch's NrnThread->data and CVODES y/ydot
  static void CopyState(Branch *b, N_Vector y, const CopyOp op);

  /// CVODES BDF max-order (NEURON=5)
  const static int kBDFMaxOrder = 5;

  /// CVODES Relative torelance (NEURON=1e-3 or 1e-4)
  constexpr static double kRelativeTolerance = 1e-3;

  /// CVODES Absolute tolerance for voltage values (NEURON=1e-8)
  constexpr static double kAbsToleranceVoltage = 1e-6;

  /// CVODES Absolute tolerance for mechanism states (NEURON=1e-8)
  constexpr static double kAbsToleranceMechStates = 1e-3;

  /// Time-window size for grouping of events to be delivered
  /// simmultaneously (0 for no grouping, Coreneuron=0.0125)
  constexpr static double kEventsDeliveryTimeWindow = 0.0125;

  /// copy CVODES y to NrnThread->data (V and m)
  inline static void ScatterY(Branch *branch, N_Vector y);

  /// copy NrnThread->data (V and m) to CVODES y
  inline static void GatherY(Branch *branch, N_Vector y);

  /// copy CVODES ydot to NrnThread->data (RHS and dm)
  inline static void ScatterYdot(Branch *branch, N_Vector ydot);

  /// copy NrnThread->data (RHS and dm) to CVODES ydot
  inline static void GatherYdot(Branch *branch, N_Vector ydot);

  /// RHS function: given y (V and m) returns ydot (RHS and dm)
  static int RHSFunction(floble_t t, N_Vector y_, N_Vector ydot,
                         void *user_data);

  /// root function: computes g_i(t,y) for detection of AP-threshold
  static int RootFunction(realtype t, N_Vector y_, realtype *gout,
                          void *user_data);

  /// jacobian function: compute J(t,y) on a dense matrix
  static int JacobianDense(long int N, floble_t t, N_Vector y_, N_Vector fy,
                           DlsMat J, void *user_data, N_Vector tmp1,
                           N_Vector tmp2, N_Vector tmp3);

  /// Solve-function for Neuron-based diagonal solver
  static int PreConditionedDiagonalSolver(CVodeMem m, N_Vector b,
                                          N_Vector weight, N_Vector ycur,
                                          N_Vector fcur);
};

};  // interpolators

};  // neurox
