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

    static void coreNeuronFinitialize();
    static void compareMechanismsFunctionPointers( std::list<NrnThreadMembList*> & uniqueMechs);

    static void registerHpxActions();  ///> Register all HPX actions
    static hpx_action_t compareBranch; ///> compares a branch to CoreNeuron data structures
    static void compareBranch2(Branch * branch); ///compares a branch to Coreneuron data structures

    //Interfaces to CoreNeuron methods
    static void finitialize();   ///> calls Coreneuron's finitialize
    static void fixed_step_minimal(); ///> calls fadvance_core.c::nrn_fixed_step_thread()
    static void fixed_step_minimal(NrnThread * nth, int secondorder); ///> calls fixed_step_minimal for single NrnThread;
    static void stepAfterStepComparison(Branch *b, NrnThread * nth, int secondorder);
    static void stepAfterStepFinitialize(Branch *b, NrnThread *nth);

  private:
    static int compareBranch_handler();
};

}; //Coreneuron
}; //Input
}; //NeuroX
