#pragma once

#include "neurox/Neurox.h"
#include "coreneuron/nrnoc/membfunc.h"
#include "coreneuron/nrnoc/membdef.h"
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrnconf.h"
#include <queue>

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

    ///nodes
    short int n;			///> number of compartments
    double * a;				///> top diagonal of Linear Algebra sparse tridiagonal matrix
    double * b;				///> bottom diagonal of Linear Algebra sparse tridiagonal matrix
    double * d;				///> main diagonal of Linear Algebra spart tridiagonal matrix
    double * v;				///> current voltage per compartment
    double * rhs;			///> right-hand side (solution vector) of Linear Algebra solver
    double * area;			///> current area per compartment

    struct MechanismInstances
    {
        double * data;			   ///> all double data used by all mechanisms
        int * dataOffsets;         ///> offset of each mechanims type in data
        int * pdata;			   ///> pointer data (offsets) for mechanisms
        int * pdataOffsets;        ///> offset of each mechanims type in Pointer-data
        int * nodesIndices;        ///> nodeindices contains all nodes this extension is responsible for, ordered according to the matrix
        int * nodesIndicesOffsets; ///> compartments' indices for each mech type
    } mechsInstances;

    //List of branches
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
    static hpx_action_t queueSpike; ///> add incoming synapse to queue
    static hpx_action_t deliverSpikes; ///> delivering of (queued) synapses (runs NET_RECEIVE on mod files)
    static hpx_action_t gaussianBackTriangulation; ///> Gaussian elimination's back triangulation: solve_core.c:triang()
    static hpx_action_t gaussianFwdSubstitution; ///> Gaussian elimination's forward substitution: solve_core.c:bksub()
    static hpx_action_t secondOrderCurrent; ///> Second Order Current : eion.c:second_order_cur()

    ///queue of incoming spikes (delivered at the end of the step)
    std::priority_queue<Synapse> synapsesQueue;

  private:

    hpx_t synapsesQueueMutex;   ///> mutex to protect the memory access to synapses queue

    static int setupMatrixRHS_handler(const char isSoma, const double v_parent);
    static int setupMatrixLHS_handler(const char isSoma);
    static int updateV_handler(const int secondOrder);
    static int gaussianBackTriangulation_handler(const char isSoma);
    static int gaussianFwdSubstitution_handler(const char isSoma, const double parentRHS);
    static int secondOrderCurrent_handler();
    static int setupMatrixInitValues_handler();
    static int setV_handler(const double v);
    static int callMechsFunction_handler(const Mechanism::Functions functionId, const double t, const double dt);
    static int queueSpike_handler(const Synapse * syn, size_t);
    static int deliverSpikes_handler();
    static int init_handler(const int n, const double *a, const double *b, const double *d,
                            const double *v, const double *rhs, const double *area,
                            const int * mechsDataOffsets, const double *data,
                            const Datum *pdata, const int branchesCount, const hpx_t * branches);
};

}; //namespace
