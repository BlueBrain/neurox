#pragma once

#include "neurox/Neurox.h"
#include "coreneuron/nrnoc/membfunc.h"
#include "coreneuron/nrnoc/membdef.h"
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrnconf.h"
#include <queue>
#include <map>

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
    Branch()=delete;
    Branch(const int n, const double *a, const double *b, const double *d,
           const double *v, const double *rhs, const double *area,
           const int * mechsCounts, const double *data,
           const int *pdata, const int branchesCount, const hpx_t * branches);
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
        int instancesCount;        ///> number of instances of particular mechanism
        double * data;			   ///> and mechanisms data (double)
        int * pdata;			   ///> pointer data (offsets)
        int * nodesIndices;        ///> nodeindices contains all nodes this extension is responsible for, ordered according to the matrix
    } * mechsInstances;            ///> Arrays of mechanism instances (per mechanism type)

    //List of children branches
    int branchesCount;		///> number of branches
    hpx_t *branches;		///> hpx address of the branches branches

    static void registerHpxActions(); ///> Register all HPX actions
    static hpx_action_t init; ///> Initializes Branch
    static hpx_action_t setupMatrixRHS; ///> finitialize.c::nrn_finitialize
    static hpx_action_t setupMatrixLHS; ///> finitialize.c::nrn_finitialize
    static hpx_action_t updateV; ///> fadvance_core.c : update()
    static hpx_action_t setupMatrixInitValues; ///> set D and RHS of all compartments to 0
    static hpx_action_t setV; ///> finitialize.c :: sets initial values of V
    static hpx_action_t callMechsFunction; ///> calls MOD functions, and BAMembList (nrn_ba)
    static hpx_action_t queueSpikes; ///> add incoming synapse to queue
    static hpx_action_t deliverNetEvents; ///> delivering of (queued) synapses (runs NET_RECEIVE on mod files)
    static hpx_action_t gaussianBackTriangulation; ///> Gaussian elimination's back triangulation: solve_core.c:triang()
    static hpx_action_t gaussianFwdSubstitution; ///> Gaussian elimination's forward substitution: solve_core.c:bksub()
    static hpx_action_t secondOrderCurrent; ///> Second Order Current : eion.c:second_order_cur()

    /** map of incoming netcons per pre-synaptic id
     *  (Equivalent to a PreSyn for local transmission and InputPreSyn
     *  for network transmission in Coreneuron (netcon.h)) */
    std::map<int, std::vector<NetConX> > netcons;

    ///queue of incoming spikes (to be delivered at the end of the step)
    std::priority_queue<Spike> spikesQueue;

  private:

    hpx_t spikesQueueMutex;   ///> mutex to protect the memory access to synapsesQueue

    static int setupMatrixRHS_handler(const char isSoma, const double v_parent);
    static int setupMatrixLHS_handler(const char isSoma);
    static int updateV_handler(const int secondOrder);
    static int gaussianBackTriangulation_handler(const char isSoma);
    static int gaussianFwdSubstitution_handler(const char isSoma, const double parentRHS);
    static int secondOrderCurrent_handler();
    static int setupMatrixInitValues_handler();
    static int setV_handler(const double v);
    static int callMechsFunction_handler(const Mechanism::Functions functionId, const double t, const double dt);
    static int queueSpikes_handler(const int preNeuronId, double deliveryTime);
    static int deliverNetEvents_handler();
    static int init_handler(const int n, const double *a, const double *b, const double *d,
                            const double *v, const double *rhs, const double *area,
                            const int * mechsDataOffsets, const double *data,
                            const Datum *pdata, const int branchesCount, const hpx_t * branches);
};

}; //namespace
