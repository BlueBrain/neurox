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
    static void createNeuron(int gid, Compartment & topCompartment, vector<Mechanism> & mechanisms, double APThreshold, vector<Synapse> & synapses);
    static hpx_t createBranch( Compartment * topCompartment, vector<Mechanism> & mechanisms);
};

}; //Input
}; //Neurox
