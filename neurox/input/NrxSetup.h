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
class NrxSetup
{
  public:
    NrxSetup();
    ~NrxSetup();

    static void copyFromCoreneuronToHpx(); ///> Copies Coreneuron data structs to HPX

  private:
    static void createBrain(int neuronsCount);
    static void createNeuron(int gid, vector<Compartment> & compartments, vector<Mechanism> & mechanisms);
    static hpx_t createBranch( Compartment * topCompartment, vector<Mechanism> & mechanisms);
};
