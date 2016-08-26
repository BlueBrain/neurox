#pragma once

#include "neurox/neurox.h"

#include <map>
#include <vector>
#include <memory>

using namespace std;

namespace neurox
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

    static void coreNeuronFakeSteps();
    static void loadData(int argc, char ** argv); ///> Copies Coreneuron data structs to HPX

    ///converts local branch data to coreneuron compatible (to use on mechanisms calls)
    static void fromHpxToCoreneuronDataStructs(const void * branch, Memb_list & membList,
                                               NrnThread & nrnThread, int mechType);

    static void compareDataStructuresWithCoreNeuron(Branch * branch);

    static void registerHpxActions();

  private:
    static void addNetConsForThisNeuron(int neuronId, int preNeuronId, int netconsCount,
                                        int netconsOffset, map< int, vector<NetConX*> > & netcons);
    static void coreNeuronInitialSetup(int argc, char ** argv);
    static hpx_t createBranch(hpx_t target, deque<Compartment*> & compartments, Compartment * topCompartment,
                              map< int, vector<NetConX*> > & netcons, int totalN, map<int, pair<int,int>> & offsetToInstance);

    static int getNeuronIdFromNrnThreadId(int nrn_id);
    static void getMechTypeAndInstanceForBranch(int & mechType, int & mechInstance);

    static int getBranchData(deque<Compartment*> & compartments, vector<double> & data, vector<int> & pdata, vector<void*> & vdata,
                             vector<int> & p, vector<int> & instancesCount, vector<int> & nodesIndices,
                             int totalN, map<int, pair<int,int>> & offsetToInstance);

    static void getVecPlayBranchData(deque<Compartment*> & compartments, vector<double> & vecPlayTdata,
                                    vector<double> & vecPlayYdata, vector<PointProcInfo> & vecPlayInfo);

    static void getNetConsBranchData(deque<Compartment*> & compartments, map<int, vector<NetConX*> > & netcons,
                                     vector<NetConX> & branchNetCons, vector<int> & branchNetConsPreNeuronId,
                                     vector<double> & branchNetConsArgs);

    static void printSubClustersToFile(FILE * fileCompartments, Compartment *topCompartment);

    static hpx_action_t createNeuron;
    static hpx_action_t broadcastNetCons;
    static hpx_action_t addSynapseTarget;

    static int createNeuron_handler(const int * i, const size_t);
    static int broadcastNetCons_handler(const int nargs, const void *args[], const size_t sizes[]);
    static int addSynapseTarget_handler(const hpx_t * synapseTarget, const size_t);
};

}; //Coreneuron
}; //Input
}; //NeuroX
