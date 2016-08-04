#pragma once

#include "neurox/Neurox.h"
#include <queue>
#include <map>
#include <vector>

namespace  Neurox {

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
    ~Branch();

    //sparse matrix information:
    short int n;	///> number of compartments
    double * a;		///> top diagonal of Linear Algebra sparse tridiagonal matrix
    double * b;		///> bottom diagonal of Linear Algebra sparse tridiagonal matrix
    double * d;		///> main diagonal of Linear Algebra spart tridiagonal matrix
    double * v;		///> current voltage per compartment
    double * rhs;	///> right-hand side (solution vector) of Linear Algebra solver
    double * area;	///> current area per compartment

    struct MechanismInstances
    {
        int instancesCount; ///> number of instances of particular mechanism
        double * data;	    ///> and mechanisms data (double)
        int * pdata;		///> pointer data (offsets) //TODO this could be unsigned short or int
        int * nodesIndices; ///> index of node this instance will be applied to //TODO this could be unsigned short
    } * mechsInstances;     ///> Arrays of mechanism instances (total size of Neuron::mechanismsCount)

    //List of children branches
    int branchesCount;		///> number of branches
    hpx_t *branches;		///> hpx address of the branches branches

    /** map of incoming netcons per pre-synaptic id
     *  (Equivalent to a PreSyn for local transmission and InputPreSyn
     *  for network transmission in Coreneuron (netcon.h)) */
    std::map<int, std::vector<NetConX> > netcons;

    ///queue of incoming spikes (to be delivered at the end of the step), sorted by delivery time
    std::priority_queue<Spike> spikesQueue;

    /// returns all instances of mechanisms type 'type'
    inline MechanismInstances & getMechanismInstanceFromType(int type);

    static void registerHpxActions(); ///> Register all HPX actions
    static hpx_action_t init; ///> Initializes the diagonal matrix and children branches for this branch
    static hpx_action_t initMechanismsInstances; ///> Initializes applications of mechanisms to this branch
    static hpx_action_t initNetCons; ///> Initializes Network Connections (NetCons) for this branch
    static hpx_action_t setupMatrixRHS; ///> finitialize.c::nrn_finitialize
    static hpx_action_t setupMatrixLHS; ///> finitialize.c::nrn_finitialize
    static hpx_action_t updateV; ///> fadvance_core.c : update()
    static hpx_action_t setupMatrixInitValues; ///> set D and RHS of all compartments to 0
    static hpx_action_t setV; ///> finitialize.c :: sets initial values of V
    static hpx_action_t callModFunction; ///> calls MOD functions, and BAMembList (nrn_ba)
    static hpx_action_t callNetReceiveFunction; ///> calls NetReceive Functions
    static hpx_action_t queueSpikes; ///> add incoming synapse to queue
    static hpx_action_t gaussianBackTriangulation; ///> Gaussian elimination's back triangulation: solve_core.c:triang()
    static hpx_action_t gaussianFwdSubstitution; ///> Gaussian elimination's forward substitution: solve_core.c:bksub()
    static hpx_action_t secondOrderCurrent; ///> Second Order Current : eion.c:second_order_cur()
    static hpx_action_t getSomaVoltage; ///>returns the voltage on the first compartment of this branch (soma if top branch)

  private:

    hpx_t spikesQueueMutex;   ///> mutex to protect the memory access to spikesQueue

    static int setupMatrixRHS_handler(const int nargs, const void *args[], const size_t sizes[]);
    static int setupMatrixLHS_handler(const char * isSoma, const size_t);
    static int updateV_handler(const int * secondOrder, const size_t);
    static int gaussianBackTriangulation_handler(const char * isSoma, const size_t);
    static int gaussianFwdSubstitution_handler(const int nargs, const void *args[], const size_t sizes[]);
    static int secondOrderCurrent_handler();
    static int setupMatrixInitValues_handler();
    static int setV_handler(const double * v, const size_t);
    static int callNetReceiveFunction_handler(const int nargs, const void *args[], const size_t sizes[]);
    static int callModFunction_handler(const Mechanism::ModFunction * functionId, const size_t);
    static int queueSpikes_handler(const int nargs, const void *args[], const size_t sizes[]);
    static int init_handler(const int nargs, const void *args[], const size_t sizes[]);
    static int initMechanismsInstances_handler(const int nargs, const void *args[], const size_t sizes[]);
    static int initNetCons_handler(const int nargs, const void *args[], const size_t sizes[]);
    static int getSomaVoltage_handler();
};

}; //namespace
