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
    static hpx_t createBranch(NrnThread * nt, map<int, list<int> > & tree, map<int, list < tuple< int, double*, int*> > > & mechanisms, int topNodeId );
    static void createNeuron(NrnThread & nt, int nrnThreadId, int neuronId);
};
