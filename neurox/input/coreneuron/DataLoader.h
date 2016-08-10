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
    static void fromHpxToCoreneuronDataStructs(void * branch, Memb_list & membList,
                                               NrnThread & nrnThread, int mechType);

    static void compareDataStructuresWithCoreNeuron(Branch * branch);
  private:
    static void addNetConsForThisNeuron(int neuronId, int preNeuronId, int netconsCount,
                                        int netconsOffset, map<int, std::vector<NetConX*> > & netcons);
    static void coreNeuronInitialSetup(int argc, char ** argv);
    static hpx_t createBranch(char isSoma, vector<Compartment*> & compartments, Compartment * topCompartment,  map<int, vector<NetConX*> > & netcons);

    static int getNeuronIdFromNrnThreadId(int nrn_id);
    static void getMechTypeAndInstanceForBranch(int & mechType, int & mechInstance);

private:
    static Compartment* getBranchingMultispliX(Compartment * topCompartment, vector<double> & d, vector<double> & b,
                                               vector<double> & a, vector<double> & rhs, vector<double> & v, vector<double> & area,
                                               vector<int> & p, vector<int> & instancesCount, vector<vector<double>> & data,
                                               vector<vector<int>> & pdata, vector<vector<int>> & nodesIndices);

    static void getBranchingFlat(vector<Compartment*> & compartments, vector<double> & data, vector<int> & pdata,
                                 vector<int> & p, vector<int> & instancesCount, vector<int> & nodesIndices);
};

}; //Coreneuron
}; //Input
}; //NeuroX
