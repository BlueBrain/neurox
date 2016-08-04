#pragma once

#include "neurox/Neurox.h"

#include <map>
#include <vector>
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
class DataLoader
{
  public:
    DataLoader()=delete;
    ~DataLoader()=delete;

    static void loadData(int argc, char ** argv); ///> Copies Coreneuron data structs to HPX

    ///converts local branch data to coreneuron compatible (to use on mechanisms calls)
    static void fromHpxToCoreneuronDataStructs(Branch * branch, Memb_list & membList,
                                               NrnThread & nrnThread, int mechType);

  private:
    static void addNetConsForThisNeuron(int neuronId, int preNeuronId, int netconsCount,
                                        int netconsOffset, map<int, std::vector<NetConX*> > & netcons);
    static void coreNeuronInitialSetup(int argc, char ** argv);
    static hpx_t createBranch( char isSoma, Compartment * topCompartment, map<int, vector<NetConX*> > & netcons);

    static int getNeuronIdFromNrnThreadId(int nrn_id);
    static void getMechTypeAndInstanceForBranch(int & mechType, int & mechInstance);
};

}; //Coreneuron
}; //Input
}; //Neurox
