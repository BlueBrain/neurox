#pragma once

#include "neurox/neurox.h"

#include <map>
#include <vector>
#include <memory>

#define OUTPUT_COMPARTMENTS_NRNTHREAD_DOT_FILE false
#define DOT_PNG_BACKGROUND_COLOR "white" //"transparent"
#define OUTPUT_NETCONS_DOT_FILE_INCLUDE_OTHERS false

using namespace std;

namespace neurox {

namespace Input {

namespace Coreneuron {

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
    static void cleanData(); ///>removes all Nrn data structures

    static void compareDataStructuresWithCoreNeuron(Branch * branch);

    static void registerHpxActions();

  private:
    static void addNetConsForThisNeuron(
            neuron_id_t neuronId, neuron_id_t preNeuronId, int netconsCount,
            offset_t netconsOffset, map< neuron_id_t, vector<NetConX*>> & netcons);

    static hpx_t createBranch(
            hpx_t target, deque<Compartment*> & compartments,
            Compartment * topCompartment,
            map< neuron_id_t, vector<NetConX*>> & netcons,
            int totalN, map<offset_t, pair<int,offset_t>> & offsetToInstance);

    static neuron_id_t getNeuronIdFromNrnThreadId(int nrn_id);
    static void getMechTypeAndInstanceForBranch(
            int & mechType, int & mechInstance);

    static int getBranchData(
            deque<Compartment*> & compartments, vector<floble_t> & data,
            vector<offset_t> & pdata, vector<void*> & vdata,
            vector<offset_t> & p, vector<offset_t> & instancesCount,
            vector<offset_t> & nodesIndices, int totalN,
            map<offset_t, pair<int,offset_t>> & offsetToInstance);

    static void getVecPlayBranchData(
            deque<Compartment*> & compartments, vector<floble_t> & vecPlayTdata,
            vector<floble_t> & vecPlayYdata, vector<PointProcInfo> & vecPlayInfo);

    static void getNetConsBranchData(
            deque<Compartment*> & compartments,
            map<neuron_id_t, vector<NetConX*>> & netcons,
            vector<NetConX> & branchNetCons,
            vector<neuron_id_t> & branchNetConsPreId,
            vector<floble_t> & branchNetConsArgs);

    static void printSubClustersToFile(
            FILE * fileCompartments, Compartment *topCompartment);

    static hpx_action_t createNeuron;
    static hpx_action_t broadcastNetCons;
    static hpx_action_t addSynapseTarget;

    static int createNeuron_handler(const int * i, const size_t);
    static int broadcastNetCons_handler(const int, const void *[],const size_t[]);
    static int addSynapseTarget_handler(const hpx_t * synapseTarget, const size_t);
};

}; //Coreneuron
}; //Input
}; //NeuroX
