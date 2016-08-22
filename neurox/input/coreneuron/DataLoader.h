#pragma once

#include "neurox/Neurox.h"

#include <map>
#include <vector>
#include <memory>

using namespace std;

namespace NeuroX
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
  private:
    static void addNetConsForThisNeuron(int neuronId, int preNeuronId, int netconsCount,
                                        int netconsOffset, vector< vector<NetConX*> > & netcons);
    static void coreNeuronInitialSetup(int argc, char ** argv);
    static hpx_t createBranch(char isSoma, deque<Compartment*> & compartments, Compartment * topCompartment,
                              vector< vector<NetConX*> > & netcons, int totalN, map<int, pair<int,int>> & offsetToInstance);

    static int getNeuronIdFromNrnThreadId(int nrn_id);
    static void getMechTypeAndInstanceForBranch(int & mechType, int & mechInstance);

private:
    static int getBranchData(deque<Compartment*> & compartments, vector<double> & data, vector<int> & pdata, vector<void*> & vdata,
                             vector<int> & p, vector<int> & instancesCount, vector<int> & nodesIndices,
                             int totalN, map<int, pair<int,int>> & offsetToInstance);

    static int getVecPlayBranchData(deque<Compartment*> & compartments, vector<double> & vecPlayTdata,
                                    vector<double> & vecPlayYdata, vector<PointProcInfo> & vecPlayInfo);

#ifdef DEBUG
    static void printSubClustersToFile(FILE * fileCompartments, Compartment *topCompartment);
#endif
};

}; //Coreneuron
}; //Input
}; //NeuroX
