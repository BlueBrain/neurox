#pragma once

#include "neurox/neurox.h"

#include <map>
#include <tuple>
#include <list>
#include <set>
#include <memory>

using namespace std;

/**
 * @brief The NrxSetup class
 * Converts CoreNeuron data structures to HPX data structures
 */
class CoreNeuronDataLoader //: IDataLoader
{
  public:
    CoreNeuronDataLoader();
    ~CoreNeuronDataLoader();

    static void loadData(); ///> Copies Coreneuron data structs to HPX

  private:
    static void createBrain(int neuronsCount, vector<Mechanism> & mechanisms);
    static void createNeuron(int gid, Compartment & topCompartment, vector<Mechanism> & mechanisms);
    static hpx_t createBranch( Compartment * topCompartment, vector<Mechanism> & mechanisms);
};
