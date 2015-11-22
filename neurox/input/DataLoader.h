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

namespace input {

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
    static void InitAndLoadCoreneuronData(
            int argc, char ** argv,
            bool nrnmpi_under_nrncontrol=false,
            bool run_setup_cleanup=false); ///> call coreneuron nrn_init_and_load_data
    static void CleanCoreneuronData(const bool clean_ion_global_map = false); ///>removes all Nrn data structures
    static void RegisterHpxActions();

    static hpx_action_t Init;
    static hpx_action_t InitMechanisms;
    static hpx_action_t InitNeurons;
    static hpx_action_t setNeurons;
    static hpx_action_t InitNetcons;
    static hpx_action_t Finalize;

    class IonInstancesInfo
    {
      public:
        int mechType;
        offset_t dataStart, dataEnd; ///> beginning and end offsets of instances in data
        vector<int> nodeIds; ///> compartments ids for each instance
    };

  private:

    static hpx_t CreateBranch(int nrnThreadId,  hpx_t somaBranchAddr,
                              deque<Compartment*> & allCompartments, Compartment * topCompartment,
                              vector<DataLoader::IonInstancesInfo> & ionsInstancesInfo, int branchingDepth,
                              int thvar_index=-1 /*AIS*/, floble_t APthreshold=0 /*AIS*/);

    static neuron_id_t GetNeuronIdFromNrnThreadId(int nrn_id);
    static void getMechTypeAndInstanceForBranch(
            int & mechType, int & mechInstance);

    static int GetBranchData(
            deque<Compartment*> & compartments, vector<floble_t> & data,
            vector<offset_t> & pdata, vector<unsigned char> & vdata,
            vector<offset_t> & p, vector<offset_t> & instancesCount,
            vector<offset_t> & nodesIndices, int N,
            vector<DataLoader::IonInstancesInfo> & ionsInstancesInfo,
            vector<map<int,int>> * mechInstanceMap = NULL);

    static void GetVecPlayBranchData(
            deque<Compartment*> & compartments, vector<floble_t> & vecPlayTdata,
            vector<floble_t> & vecPlayYdata, vector<PointProcInfo> & vecPlayInfo,
            vector<map<int,int>> * mechInstanceMap = NULL);

    static void GetNetConsBranchData(
            deque<Compartment*> & compartments,
            vector<NetConX> & branchNetCons,
            vector<neuron_id_t> & branchNetConsPreId,
            vector<floble_t> & branchNetConsArgs,
            vector<map<int,int>> * mechInstanceMap = NULL);

    static void GetAllChildrenCompartments(
            deque<Compartment*> & subSection,
            Compartment * topCompartment);

    static void GetMechInstanceMap(
            deque<Compartment*> & compartments,
            vector<map<int,int>> & mechsInstancesMap);


    static void SetMechanisms2(
            const int mechsCount, const int* mechsIds,
            const int* dependenciesCount, const int* dependencies,
            const int* successorsCount, const int* successors); ///> Set Mechanisms

    static void PrintSubClustersToFile(
            FILE * fileCompartments, Compartment *topCompartment);

    static PointProcInfo GetPointProcInfoFromDataPointer(NrnThread * nt, double *pd, size_t size);

    static hpx_action_t AddSynapse;
    static hpx_action_t AddNeurons;
    static hpx_action_t SetMechanisms;

    static int CreateNeuron(int neuron_idx, void *);
    static int GetMyNrnThreadsCount();
    static int AddSynapse_handler(const int, const void *[], const size_t[]) ;
    static int AddNeurons_handler(const int, const void *[], const size_t[]) ;
    static int SetMechanisms_handler(const int, const void *[], const size_t[]) ;
    static int Init_handler ();
    static int InitMechanisms_handler();
    static int InitNeurons_handler();
    static int setNeurons_handler();
    static int InitNetcons_handler();
    static int Finalize_handler();
};

}; //Input
}; //NeuroX

