#pragma once

#include "neurox/neurox.h"

#include <list>
#include <memory>
#include <vector>

using namespace std;
namespace neurox {
namespace input {

/**
 * @brief The DataComparison class
 * Compares CoreNeuron and HPX-based data structures
 */
class Debugger {
 public:
  Debugger() = delete;
  ~Debugger() = delete;

  static bool IsEqual(floble_t a, floble_t b, bool roughlyEqual = false);

  static void CompareAllBranches();  ///> compares all branches wth Coreneuron
  static void CompareBranch2(
      Branch *branch);  ///> compares a branch to Coreneuron data structures
  static void CompareMechanismsFunctions();

  /// compares a branch to CoreNeuron data structures
  static hpx_action_t CompareBranch;

  /// calls finitialize (on a compute node, not Branch)
  static hpx_action_t Finitialize;

  /// calls nrn_fixed_step_group_minimal (on a compute node, not Branch)
  static hpx_action_t NrnSpikeExchange;

  static hpx_action_t FixedStepMinimal;
  static hpx_action_t ThreadTableCheck;

  // Interfaces to CoreNeuron methods

  /// calls fixed_step_minimal for single NrnThread;
  static void FixedStepMinimal2(NrnThread *nth, int secondorder);
  static void StepAfterStepBackwardEuler(Branch *b, NrnThread *nth,
                                         int secondorder);
  static void StepAfterStepFinitialize(Branch *b, NrnThread *nth);
  static void RunCoreneuronAndCompareAllBranches();
  static void SingleNeuronStepAndCompare(NrnThread *nt, Branch *b,
                                         char secondorder);

  static void RegisterHpxActions();  ///> Register all HPX actions

 private:
  static int CompareBranch_handler();
  static int Finitialize_handler();
  static int NrnSpikeExchange_handler();
  static int FixedStepMinimal_handler(const int *, const size_t);
  static int ThreadTableCheck_handler();
};
};  // namespace input
};  // namespace neurox
