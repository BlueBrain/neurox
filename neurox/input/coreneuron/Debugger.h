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

    static void coreNeuronFinitialize();
    static void compareMechanismsFunctionPointers( std::list<NrnThreadMembList*> & uniqueMechs);

    static void registerHpxActions();  ///> Register all HPX actions
    static hpx_action_t compareBranch; ///> compares a branch to CoreNeuron data structures
    static void finitialize();   ///> calls Coreneuron's finitialize

  private:
    static int compareBranch_handler();
};

}; //Coreneuron
}; //Input
}; //NeuroX
