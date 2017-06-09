#pragma once

#include "neurox/neurox.h"

#include <list>
#include <vector>
#include <memory>

using namespace std;
namespace neurox {
namespace Input {
namespace Coreneuron {

/**
 * @brief The DataComparison class
 * Compares CoreNeuron and HPX-based data structures
 */
class Debugger
{
  public:
    Debugger()=delete;
    ~Debugger()=delete;

    static bool isEqual(floble_t a, floble_t b, bool roughlyEqual = false);

    static void compareAllBranches(); ///> compares all branches wth Coreneuron
    static void compareBranch2(Branch * branch); ///> compares a branch to Coreneuron data structures
    static void compareMechanismsFunctionPointers();

    static hpx_action_t compareBranch; ///> compares a branch to CoreNeuron data structures
    static hpx_action_t finitialize;   ///> calls finitialize (on a compute node, not Branch)
    static hpx_action_t nrnSpikeExchange;   ///> calls nrn_fixed_step_group_minimal (on a compute node, not Branch)
    static hpx_action_t fixedStepMinimal;
    static hpx_action_t threadTableCheck;

    //Interfaces to CoreNeuron methods
    static void fixed_step_minimal(NrnThread * nth, int secondorder); ///> calls fixed_step_minimal for single NrnThread;
    static void stepAfterStepBackwardEuler(Branch *b, NrnThread * nth, int secondorder);
    static void stepAfterStepFinitialize(Branch *b, NrnThread *nth);

    static void registerHpxActions();  ///> Register all HPX actions

  private:
    static int compareBranch_handler();
    static int finitialize_handler();
    static int nrnSpikeExchange_handler();
    static int fixedStepMinimal_handler(const int *, const size_t);
    static int threadTableCheck_handler();
};

}; //Coreneuron
}; //Input
}; //NeuroX
