#pragma once

#include "neurox/neurox.h"
#include <queue>
#include <map>
#include <vector>
#include <set>
#include <new>  //placement new

using namespace neurox;

namespace  neurox {

class Neuron;

/**
 * @brief The Branch class
 * Represents a branch as a continuous set of compartments (nodes)
 * and a bifurcation made of branches branchs spreading out from the last node;
 * Also represents mechanisms as a list of mechanisms and their applications
 * to several nodes of this branch.
 */
class Branch
{
  public:
    Branch()=delete; ///> no constructor, build using hpx init function instead
    static void* operator new(size_t bytes, void* addr);
    static void  operator delete(void* worker);

    Branch(offset_t n,
           hpx_t branchHpxAddr,
           floble_t * data, size_t dataCount,
           offset_t *pdata, size_t pdataCount,
           offset_t * instancesCount, size_t instancesCountCount, //same tipe as pdata
           offset_t * nodesIndices, size_t nodesIndicesCount,
           hpx_t * branches, size_t branchesCount,
           offset_t * p, size_t pCount,
           floble_t * vecplayT, size_t vecplayTCount,
           floble_t * vecplayY, size_t vecplayYCount,
           PointProcInfo * vecplayInfo, size_t vecplayCount,
           NetConX * netcons, size_t netconsCount,
           neuron_id_t * netConsPreId, size_t netConsPreIdsCount,
           floble_t *netConsArgs, size_t netConsArgsCount,
           void** vdata, size_t vdataCount);
    ~Branch();

    NrnThread *nt;             ///> compartments metadata
    Memb_list *mechsInstances; ///> Arrays of mechanism instances (size neurox::mechanismsCount)
    Neuron * soma;             ///> if top branch, it's populated, otherwise is NULL

    struct MechanismsGraphLCO
    {
        hpx_t * mechsLCOs; ///>contains the HPX address of the and-gate of each mechanism in the graph
        hpx_t endLCO; ///> represents the bottom of the graph
        hpx_t graphLCO; ///> controls all active thread on the graph of mechanisms (including the 'end' node)

        static hpx_action_t nodeFunction; ///> represents the action of the nodes in the mechanisms graph
        static int nodeFunction_handler(const int * mechType_ptr, const size_t);
    } * mechsGraph; ///> represents the parallel computation graph of mechanisms instances (NULL for serial)

    void initMechanismsGraph(hpx_t target); ///> creates mechanisms instance graph based on global var 'mechanisms'

    struct NeuronTree
    {
        size_t branchesCount;	///> number of branches (>0)
        hpx_t *branches;		///> hpx address of the branches branches

        static const size_t futuresSize = 3; ///> size of futures arrays (used in Gaussian elimination)
        hpx_t localLCO[futuresSize]; ///> LCO of current branch execution to communicate 3 variables with parent branch
        hpx_t (*branchesLCOs)[futuresSize]; ///> LCOs of branches' executiont (NULL if no children)
    } * neuronTree; ///> represents the tree structure (or NULL if none i.e. full neuron representation)

    std::map<neuron_id_t, std::vector<NetConX*> > netcons; ///> map of incoming netcons per pre-synaptic id

    std::priority_queue< std::pair<floble_t,Event*> > eventsQueue;  ///>queue of incoming events sorted per delivery time
    hpx_t eventsQueueMutex;   ///> mutex to protect the memory access to eventsQueue

    static void registerHpxActions(); ///> Register all HPX actions
    static hpx_action_t init; ///> Initializes the diagonal matrix and children branches for this branch
    static hpx_action_t initSoma; ///> Initializes soma information in this branch
    static hpx_action_t initNeuronTreeLCO; ///> Initializes neuronTree
    static hpx_action_t clear; ///> deletes all data structures in branch and sub-branches
    static hpx_action_t addSpikeEvent; ///> add incoming synapse to queue
    static hpx_action_t finitialize;
    static hpx_action_t backwardEuler;
    static hpx_action_t backwardEulerStep;

    void callModFunction(const Mechanism::ModFunction functionId);
    void initEventsQueue(); ///> start NetEvents and PlayVect on events queue
    void deliverEvents(floble_t t);
    void fixedPlayContinuous();
    void setupTreeMatrixMinimal();
    void finitialize2();
    void backwardEulerStep2();

  private:
    static int init_handler(const int, const void *[], const size_t[]);
    static int initSoma_handler(const int, const void *[], const size_t[]);
    static int initNeuronTreeLCO_handler();
    static int clear_handler();
    static int addSpikeEvent_handler(const int, const void *[], const size_t[]);
    static int finitialize_handler();
    static int backwardEuler_handler();
    static int backwardEulerStep_handler(const int*, const size_t);
};

}; //namespace
