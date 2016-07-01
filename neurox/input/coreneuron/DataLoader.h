#pragma once

#include "neurox/Neurox.h"

#include <map>
#include <tuple>
#include <list>
#include <set>
#include <memory>

using namespace std;

namespace Neurox
{

namespace Input
{

namespace Coreneuron
{

/**
 * @brief The NrxSetup class
 * Converts CoreNeuron data structures to HPX data structures
 */
class CoreNeuronDataLoader
{
  public:
    CoreNeuronDataLoader()=delete;
    ~CoreNeuronDataLoader()=delete;

    static void loadData(int argc, char ** argv); ///> Copies Coreneuron data structs to HPX

  private:
    static void coreNeuronInitialSetup(int argc, char ** argv);
    static void createBrain(int neuronsCount, Mechanism * mechanisms, int mechanismsCount);
    static void createNeuron(int gid, Compartment & topCompartment, APthreshold);
    static hpx_t createBranch( Compartment * topCompartment);

    static void getNeuronIdFromNrnThreadId(int nrn_id);
    static void getMechTypeAndInstanceForBranch(int & mechType, int & mechInstance);
};

}; //Coreneuron
}; //Input
}; //Neurox
