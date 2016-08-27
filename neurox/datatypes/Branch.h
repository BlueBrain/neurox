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
    static void operator delete(void* worker);

    Branch(int n,
           hpx_t branchHpxAddr,
           double * data, size_t dataCount,
           int *pdata, size_t pdataCount,
           int * instancesCount, size_t instancesCountCount,
           int * nodesIndices, size_t nodesIndicesCount,
           hpx_t * branches, size_t branchesCount,
           int * p, size_t pCount,
           double * vecplayT, size_t vecplayTCount,
           double * vecplayY, size_t vecplayYCount,
           PointProcInfo * vecplayInfo, size_t vecplayCount,
           NetConX * netcons, size_t netconsCount,
           int * netConsPreId, size_t netConsPreIdsCount,
           double *netConsArgs, size_t netConsArgsCount,
           void** vdata, size_t vdataCount);
    ~Branch();

    //for conveniency we add t and dt to all branch localities
    double t;       ///> current simulation time
    double dt;      ///> time step
    char secondOrder;

    Neuron * soma;  ///> if top branch, it's populated, otherwise is NULL

    //sparse matrix information:
    int n;		    ///> number of compartments
    double * a;		///> top diagonal of Linear Algebra sparse tridiagonal matrix
    double * b;		///> bottom diagonal of Linear Algebra sparse tridiagonal matrix
    double * d;		///> main diagonal of Linear Algebra spart tridiagonal matrix
    double * v;		///> current voltage per compartment
    double * rhs;	///> right-hand side (solution vector) of Linear Algebra solver
    double * area;	///> current area per compartment
    int *p;         ///> index of parents compartments (if multiSpliX is 0) or NULL (if multiSpliX is 1)

    struct MechanismInstance
    {
        int count;          ///> number of instances of particular mechanism
        double * data;	    ///> pointer to Branch::data vector with start position of this mechanism's data
        int * pdata;		///> pointer to Branch::pdata vector with start position of this mechanism's pointer data
        int * nodesIndices; ///> array of nodes this instance will be applied to
    } * mechsInstances;     ///> Arrays of mechanism instances (total size of Neuron::mechanismsCount)

    struct MechanismsGraphLCO
    {
        hpx_t * mechsLCOs; ///>contains the HPX address of the and-gate of each mechanism in the graph
        hpx_t endLCO; ///> represents the bottom of the graph
        hpx_t graphLCO; ///> controls all active thread on the graph of mechanisms (including the 'end' node)

        static hpx_action_t nodeFunction; ///> represents the action of the nodes in the mechanisms graph
        static int nodeFunction_handler(const int * mechType_ptr, const size_t);
    } * mechsGraph; ///> represents the parallel computation graph of mechanisms instances (NULL for serial)
    void initMechanismsGraph(hpx_t target); ///> creates mechanisms instance graph based on global var 'mechanisms'

    struct NeuronTreeLCO
    {
        int branchesCount;		///> number of branches (>0)
        hpx_t *branches;		///> hpx address of the branches branches

        hpx_t parentLCO;        ///> LCO of parent execution
        hpx_t localLCO;         ///> LCO of current branch execution
        hpx_t * branchesLCOs;   ///> LCOs of branches' execution
    } * neuronTree; ///> represents the tree structure (or NULL if full neuron representation)

    std::map<int, std::vector<NetConX*> > netcons; ///> map of incoming netcons per pre-synaptic id

    std::priority_queue< std::pair<double,Event*> > eventsQueue;  ///>queue of incoming events sorted per delivery time
    hpx_t eventsQueueMutex;   ///> mutex to protect the memory access to eventsQueue

    static void registerHpxActions(); ///> Register all HPX actions
    static hpx_action_t init; ///> Initializes the diagonal matrix and children branches for this branch
    static hpx_action_t initSoma; ///> Initializes soma information in this branch
    static hpx_action_t initNeuronTreeLCO; ///> Initializes neuronTree
    static hpx_action_t clear; ///> deletes all data structures in branch and sub-branches
    static hpx_action_t addSpikeEvent; ///> add incoming synapse to queue
    static hpx_action_t finitialize;
    static hpx_action_t backwardEuler;

    double * data; ///> all double data for the branch (RHS, D, A, B, V, Area, and mechanisms)
    void ** vdata; ///> TODO make part of mechsInstances. all pointer data for the branch (RNG + something for ProbAMBA and ProbGABA)
    VecPlayContinuouX ** vecplay; ///> described continuous events
    int vecplayCount; //number of vecplay

    void callModFunction(const Mechanism::ModFunction functionId);
    void initEventsQueue(); ///> start NetEvents and PlayVect on events queue
    void deliverEvents(double t);
    void fixedPlayContinuous();
    void setupTreeMatrixMinimal();
    void secondOrderCurrent();
    void initialize();
    void backwardEulerStep();

  private:
    static int init_handler(const int nargs, const void *args[], const size_t sizes[]);
    static int initSoma_handler(const int nargs, const void *args[], const size_t sizes[]);
    static int initNeuronTreeLCO_handler(const hpx_t * parentLCO_ptr, size_t);
    static int clear_handler();
    static int addSpikeEvent_handler(const int nargs, const void *args[], const size_t sizes[]);
    static int finitialize_handler();
    static int backwardEuler_handler();
};

}; //namespace
