#pragma once

#include "neurox/Neurox.h"
#include <queue>
#include <map>
#include <vector>

using namespace std;

namespace  NeuroX {

/**
 * @brief The Branch class
 * Represents a branch as a continuous set of compartments (nodes)
 * and a bifurcation made of branches branchs spreading out from the last node;
 * Also represents mechanisms as a list of mechanisms and their applications
 * to several nodes of this branch.
 */
class Branch
{
  //TODO variables n, p* and instance.count and indices can be short int unless we merge several neurons (a la CoreNeuron)

  public:
    Branch()=delete; ///> no constructor, build using hpx init function instead
    ~Branch();

    //for conveniency we add t and dt to all branch localities
    double t;       ///> current simulation time
    double dt;      ///> time step

    //sparse matrix information:
    int n;		    ///> number of compartments
    char isSoma;    ///> is this branch the top branch of the morphology tree (ie constains soma?)
    double * a;		///> top diagonal of Linear Algebra sparse tridiagonal matrix
    double * b;		///> bottom diagonal of Linear Algebra sparse tridiagonal matrix
    double * d;		///> main diagonal of Linear Algebra spart tridiagonal matrix
    double * v;		///> current voltage per compartment
    double * rhs;	///> right-hand side (solution vector) of Linear Algebra solver
    double * area;	///> current area per compartment
    int *p;         ///> index of parents compartments (if multiSpliX is 0) or NULL (if multiSpliX is 1)

    //List of children branches
    int branchesCount;		///> number of branches (if any)
    hpx_t *branches;		///> hpx address of the branches branches

    struct MechanismInstance
    {
        int count;          ///> number of instances of particular mechanism
        double * data;	    ///> pointer to Branch::data vector with start position of this mechanism's data
        int * pdata;		///> pointer to Branch::pdata vector with start position of this mechanism's pointer data
        int * nodesIndices; ///> array of nodes this instance will be applied to
    } * mechsInstances;     ///> Arrays of mechanism instances (total size of Neuron::mechanismsCount)

    struct MechanismsExecutionGraph
    {
        hpx_t * mechsLCOs; ///>contains the HPX address of the and-gate of each mechanism in the graph
        hpx_t endLCO; ///> represents the bottom of the graph
        hpx_t graphLCO; ///> controls all active thread on the graph of mechanisms (including the 'end' node)

        static hpx_action_t nodeFunction; ///> represents the action of the nodes in the mechanisms graph
        static int nodeFunction_handler(const int * mechType_ptr, const size_t);
    } * mechsGraph; ///> represents the parallel computation graph of mechanisms instances (NULL for serial)

    map<int, deque<NetConX> > netcons; ///> map of incoming netcons per pre-synaptic id

    priority_queue< pair<double,Event*> > eventsQueue;  ///>queue of incoming events sorted per delivery time
    hpx_t eventsQueueMutex;   ///> mutex to protect the memory access to spikesQueue

    static void registerHpxActions(); ///> Register all HPX actions
    static hpx_action_t init; ///> Initializes the diagonal matrix and children branches for this branch
    static hpx_action_t clear; ///> deletes all data structures in branch and sub-branches
    static hpx_action_t initNetCons; ///> Initializes Network Connections (NetCons) for this branch
    static hpx_action_t updateV; ///> fadvance_core.c : update()
    static hpx_action_t callModFunction; ///> calls MOD functions, and BAMembList (nrn_ba)
    static hpx_action_t queueSpikes; ///> add incoming synapse to queue
    static hpx_action_t secondOrderCurrent; ///> Second Order Current : eion.c:second_order_cur()
    static hpx_action_t getSomaVoltage; ///>returns the voltage on the first compartment of this branch (soma if top branch)

    static hpx_action_t fixedPlayContinuous; ///> nrn_fixed_play_continuous
    static hpx_action_t deliverEvents;   ///> delivers all events for the next time step (nrn_deliver_events)

    double * data; ///> all double data for the branch (RHS, D, A, B, V, Area, and mechanisms)
    void ** vdata; ///> TODO make part of mechsInstances. all pointer data for the branch (RNG + something for ProbAMBA and ProbGABA)
    VecPlayContinuouX ** vecplay; ///> described continuous events
    int vecplayCount; //number of vecplay

    void callModFunction2(const Mechanism::ModFunction functionId);
    void initEventsQueue(); ///> start NetEvents and PlayVect on events queue


    static int deliverEvents_handler(const double *t, const size_t);

  private:
    static int callModFunction_handler(const Mechanism::ModFunction * functionId, const size_t);
    static int updateV_handler(const int * secondOrder, const size_t);
    static int secondOrderCurrent_handler();
    static int queueSpikes_handler(const int nargs, const void *args[], const size_t sizes[]);
    static int init_handler(const int nargs, const void *args[], const size_t sizes[]);
    static int initNetCons_handler(const int nargs, const void *args[], const size_t sizes[]);
    static int getSomaVoltage_handler();
    static int clear_handler();
    static int fixedPlayContinuous_handler();
};

}; //namespace
