#pragma once

#include "neurox/neurox.h"

#include <map>
#include <tuple>
#include <list>
#include <set>
#include <memory>

using namespace std;

//TODO should be part of Neurox::Input::CoreNeuron namespace

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
    static void createBrain(int neuronsCount, vector<Mechanism> & mechanisms);
    static void createNeuron(int gid, Compartment & topCompartment, vector<Mechanism> & mechanisms, double APThreshold, vector<Synapse> & synapses);
    static hpx_t createBranch( Compartment * topCompartment, vector<Mechanism> & mechanisms);
};
