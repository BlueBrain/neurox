#pragma once

#include "neurox/neurox.h"
#include <queue>
#include <map>
#include <vector>
#include <functional> //std::greater_equal
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
           int nrnThreadId,
           int thresholdVoffset,
           hpx_t branchHpxAddr,
           floble_t * data, size_t dataCount,
           offset_t *pdata, size_t pdataCount,
           offset_t * instancesCount, size_t instancesCountCount, //same tipe as pdata
           offset_t * nodesIndices, size_t nodesIndicesCount,
           hpx_t topBranchAddr,
           hpx_t * branches, size_t branchesCount,
           offset_t * p, size_t pCount,
           floble_t * vecplayT, size_t vecplayTCount,
           floble_t * vecplayY, size_t vecplayYCount,
           PointProcInfo * vecplayInfo, size_t vecplayCount,
           NetConX * netcons, size_t netconsCount,
           neuron_id_t * netConsPreId, size_t netConsPreIdsCount,
           floble_t *branchWeights, size_t branchWeightsCount,
           unsigned char* vdataSerialized, size_t vdataSerializedCount);
    ~Branch();

    NrnThread *nt;             ///> compartments metadata
    Memb_list *mechsInstances; ///> Arrays of mechanism instances (size neurox::mechanismsCount)
    Neuron * soma;             ///> if top branch, it's populated, otherwise is NULL
    floble_t* thvar_ptr;       ///> pointer to variable holding threashold var (if any; typically RHS[0] is AIS's first compartment)

    class MechanismsGraph
    {
      public:
        MechanismsGraph(int); ///> creates mechanisms instance graph based on global var 'mechanisms'
        ~MechanismsGraph();
        void InitMechsGraph(hpx_t branchHpxAddr); ///> Launch HPX-threads for dorment mechs-graph

        hpx_t * mechsLCOs; ///>contains the HPX address of the and-gate of each mechanism in the graph
        hpx_t endLCO; ///> represents the bottom of the graph
        hpx_t graphLCO; ///> controls all active thread on the graph of mechanisms (including the 'end' node)

        static hpx_action_t Init; ///> init function for hpx_reduce of mechanisms graphs
        static hpx_action_t Reduce; ///> reduce opreationf for hpx_creduce of mechanisms graphs
        static hpx_action_t MechFunction; ///> represents the action of the nodes in the mechanisms graph

        static int MechFunction_handler(const int * mechType_ptr, const size_t);
        static void Init_handler(Mechanism::ModFunction * func_ptr, const size_t);
        static void Reduce_handler(Mechanism::ModFunction * lhs,
                                              const Mechanism::ModFunction *rhs, const size_t);

        //for current function accumulation of shadow arrays
        hpx_t rhs_d_mutex;
        hpx_t i_didv_mutex[Mechanism::Ion::size_writeable_ions];
        static void AccumulateRHSandD  (NrnThread* nt, Memb_list * ml, int, void* args);
        static void AccumulateIandDIDV (NrnThread* nt, Memb_list * ml, int, void* args);

    } * mechsGraph; ///> represents the parallel computation graph of mechanisms instances (NULL for serial)

    class BranchTree
    {
      public:
        BranchTree()=delete;
        BranchTree(hpx_t topBranchAddr, hpx_t* branches, size_t branchesCount);
        ~BranchTree();

        hpx_t topBranchAddr;    ///> hpx address of the some branch
        hpx_t *branches;		///> hpx address of children branches
        size_t branchesCount;	///> number of branches (>0)

        static constexpr size_t futuresSize = 7; ///> size of futures arrays (used in Gaussian elimination and AP threshold)
        hpx_t withParentLCO[futuresSize]; ///> LCO of current branch execution to communicate variables with parent branch
        hpx_t (*withChildrenLCOs)[futuresSize]; ///> LCOs to communication with branches execution (NULL if no children)

        static hpx_action_t InitLCOs; ///> Initializes neuronTree
        static int InitLCOs_handler();
    } * branchTree; ///> represents the tree structure (or NULL if none i.e. full neuron representation)

    std::map<neuron_id_t, std::vector<NetConX*> > netcons; ///> map of incoming netcons per pre-synaptic gid

    /// priority queue of incoming events sorted per delivery time
#ifdef USE_TIM_SPTQ
    tools::sptq_queue <TimedEvent, std::greater_equal<TimedEvent> > eventsQueue;
 #else
    std::priority_queue< TimedEvent, std::vector<TimedEvent>, std::greater_equal<TimedEvent> > eventsQueue;
#endif
    hpx_t eventsQueueMutex;   ///> mutex to protect the memory access to eventsQueue

    static hpx_action_t Init; ///> Initializes the diagonal matrix and children branches for this branch
    static hpx_action_t InitSoma; ///> Initializes soma information in this branch
    static hpx_action_t Clear; ///> deletes all data structures in branch and sub-branches
    static hpx_action_t AddSpikeEvent; ///> add incoming synapse to queue
    static hpx_action_t UpdateTimeDependency; ///>update maximum time allowed based on received dependency info
    static hpx_action_t Finitialize; ///> finitialize.c::finitialize()
    static hpx_action_t BackwardEuler;
    static hpx_action_t BackwardEulerOnLocality;
    static hpx_action_t ThreadTableCheck;

    void CallModFunction(const Mechanism::ModFunction functionId);
    void InitVecPlayContinous(); ///> start NetEvents and PlayVect on events queue
    void AddEventToQueue(floble_t t, Event * e);
    void DeliverEvents(floble_t t);
    void FixedPlayContinuous();
    void SetupTreeMatrix();
    void SolveTreeMatrix();
    void Finitialize2();
    void BackwardEulerStep();

    static void RegisterHpxActions(); ///> Register all HPX actions

  private:
    static int Init_handler(const int, const void *[], const size_t[]);
    static int InitSoma_handler(const int, const void *[], const size_t[]);
    static int Clear_handler();
    static int AddSpikeEvent_handler(const int, const void *[], const size_t[]);
    static int UpdateTimeDependency_handler(const int, const void *[], const size_t[]);
    static int Finitialize_handler();
    static int BackwardEuler_handler(const int*, const size_t);
    static int BackwardEulerOnLocality_handler(const int*, const size_t);
    static int ThreadTableCheck_handler();
};

}; //namespace
