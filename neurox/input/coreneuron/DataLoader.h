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

  private:

    static hpx_t createBranch( int nrnThreadId,
            hpx_t target, deque<Compartment*> & compartments,
            Compartment * topCompartment,
            int totalN, map<offset_t, pair<int,offset_t>> & offsetToInstance);

    static neuron_id_t getNeuronIdFromNrnThreadId(int nrn_id);
    static void getMechTypeAndInstanceForBranch(
            int & mechType, int & mechInstance);

    static int getBranchData(
            deque<Compartment*> & compartments, vector<floble_t> & data,
            vector<offset_t> & pdata, vector<unsigned char> & vdata,
            vector<offset_t> & p, vector<offset_t> & instancesCount,
            vector<offset_t> & nodesIndices, int totalN,
            map<offset_t, pair<int,offset_t>> & offsetToInstance);

    static void getVecPlayBranchData(
            deque<Compartment*> & compartments, vector<floble_t> & vecPlayTdata,
            vector<floble_t> & vecPlayYdata, vector<PointProcInfo> & vecPlayInfo);

    static void getNetConsBranchData(
            deque<Compartment*> & compartments,
            vector<NetConX> & branchNetCons,
            vector<neuron_id_t> & branchNetConsPreId,
            vector<floble_t> & branchNetConsArgs);

    static void printSubClustersToFile(
            FILE * fileCompartments, Compartment *topCompartment);

    static PointProcInfo getPointProcInfoFromDataPointer(NrnThread * nt, double *pd);

    static hpx_action_t addSynapse;
    static hpx_action_t addNeurons;

    static int createNeuron(int neuron_idx, void * targets);
    static int getMyNrnNeuronsCount();
    static int addSynapse_handler(const int, const void *[], const size_t[]) ;
    static int addNeurons_handler(const int, const void *[], const size_t[]) ;
    static int init_handler ();
    static int initMechanisms_handler();
    static int initNeurons_handler();
    static int initNetcons_handler ();
    static int finalize_handler();
};

}; //Coreneuron
}; //Input
}; //NeuroX

