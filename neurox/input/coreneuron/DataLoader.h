#pragma once

#include "neurox/neurox.h"

#include <map>
#include <vector>
#include <memory>
#include <deque>

#define DOT_PNG_BACKGROUND_COLOR "transparent"
#define NETCONS_DOT_OUTPUT_NETCONS_FROM_EXTERNAL_NEURONS false
#define COMPARTMENTS_DOT_OUTPUT_CORENEURON_STRUCTURE false

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

    static void loadData(int argc, char ** argv); ///> Copies Coreneuron data structs to HPX
    static void initAndLoadCoreneuronData(int argc, char ** argv); ///> call coreneuron nrn_init_and_load_data
    static void cleanCoreneuronData(); ///>removes all Nrn data structures
    static void registerHpxActions();

    static hpx_action_t init;
    static hpx_action_t initMechanisms;
    static hpx_action_t initNeurons;
    static hpx_action_t initNetcons;
    static hpx_action_t finalize;

    class IonInstancesInfo
    {
      public:
        int mechType;
        offset_t dataStart, dataEnd; ///> beginning and end offsets of instances in data
        vector<int> nodeIds; ///> compartments ids for each instance
    };

  private:

    static hpx_t createBranch( int nrnThreadId,
            hpx_t target, deque<Compartment*> & allCompartments,
            Compartment * topCompartment,
            vector<DataLoader::IonInstancesInfo> & ionsInstancesInfo,
            int branchingDepth=0);

    static neuron_id_t getNeuronIdFromNrnThreadId(int nrn_id);
    static void getMechTypeAndInstanceForBranch(
            int & mechType, int & mechInstance);

    static int getBranchData(
            deque<Compartment*> & compartments, vector<floble_t> & data,
            vector<offset_t> & pdata, vector<unsigned char> & vdata,
            vector<offset_t> & p, vector<offset_t> & instancesCount,
            vector<offset_t> & nodesIndices, int totalN,
            vector<DataLoader::IonInstancesInfo> & ionsInstancesInfo,
            vector<map<int,int>> * mechInstanceMap = NULL);

    static void getVecPlayBranchData(
            deque<Compartment*> & compartments, vector<floble_t> & vecPlayTdata,
            vector<floble_t> & vecPlayYdata, vector<PointProcInfo> & vecPlayInfo,
            vector<map<int,int>> * mechInstanceMap = NULL);

    static void getNetConsBranchData(
            deque<Compartment*> & compartments,
            vector<NetConX> & branchNetCons,
            vector<neuron_id_t> & branchNetConsPreId,
            vector<floble_t> & branchNetConsArgs,
            vector<map<int,int>> * mechInstanceMap = NULL);

    static void getAllChildrenCompartments(
            deque<Compartment*> & subSection,
            Compartment * topCompartment);

    static void getMechInstanceMap(
            deque<Compartment*> & compartments,
            vector<map<int,int>> & mechsInstancesMap);

    static void printSubClustersToFile(
            FILE * fileCompartments, Compartment *topCompartment);

    static PointProcInfo getPointProcInfoFromDataPointer(NrnThread * nt, double *pd, size_t size);

    static hpx_action_t addSynapse;
    static hpx_action_t addNeurons;

    static int createNeuron(int neuron_idx, void * targets);
    static int getMyNrnNeuronsCount();
    static int addSynapse_handler(const int, const void *[], const size_t[]) ;
    static int addNeurons_handler(const int, const void *[], const size_t[]) ;
    static int init_handler ();
    static int initMechanisms_handler();
    static int initNeurons_handler();
    static int initNetcons_handler(const hpx_t* = nullptr, const size_t = 0);
    static int finalize_handler();
};

}; //Coreneuron
}; //Input
}; //NeuroX

