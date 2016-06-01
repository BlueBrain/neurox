#pragma once

#include "neurox/Neurox.h"
#include "coreneuron/nrnoc/membfunc.h"
#include "coreneuron/nrnoc/membdef.h"
#include "coreneuron/nrnconf.h"
#include <queue>

namespace  Neurox {

class FwSubFutureData;

/**
 * @brief The Branch class
 * Represents a branch as a continuous set of compartments (nodes)
 * and a bifurcation made of children branchs spreading out from the last node;
 * Also represents mechanisms as a list of mechanisms and their applications
 * to several nodes of this branch.
 */
class Branch
{
  public:
    Branch()=delete;
    Branch(const int n, const double *a, const double *b, const double *d,
           const double *v, const double *rhs, const double *area,
           const int m, const int * mechsCount, const double *data,
           const Datum *pdata, const int childrenCount, const hpx_t * children);
    ~Branch();

    //nodes
    short int n;			///> number of compartments
    double * a;				///> top diagonal of Linear Algebra sparse tridiagonal matrix
    double * b;				///> bottom diagonal of Linear Algebra sparse tridiagonal matrix
    double * d;				///> main diagonal of Linear Algebra spart tridiagonal matrix
    double * v;				///> current voltage per compartment
    double * rhs;			///> right-hand side (solution vector) of Linear Algebra solver
    double * area;			///> current area per compartment

    //Mechanisms data
    int m;                  ///> number of mechanisms applicatins
    int * mechsOffsets;     ///> offset of each mechanims type
    double * data;			///> all double data used by all mechanisms
    int * pdata;			///> all int data used by all mechanisms
    //int * pntProcTargetId;	///> when using point process, the index of the compartment it updates

    //List of children
    int childrenCount;		///> number of children branches (always >1)
    hpx_t * children;		///> hpx address of the children branches

    static void registerHpxActions(); ///> Register all HPX actions
    static hpx_action_t init; ///> Initializes Branch
    static hpx_action_t setupMatrixRHS; ///> finitialize.c::nrn_finitialize
    static hpx_action_t setupMatrixLHS; ///> finitialize.c::nrn_finitialize
    static hpx_action_t setupMatrixInitValues; ///> set D and RHS of all compartments to 0
    static hpx_action_t setV; ///> finitialize.c :: sets initial values of V
    static hpx_action_t callMechsFunction(const int functionId); ///> BAMembList and nrn_ba
    static hpx_action_t queueSpike; ///> add incoming synapse to queue
    static hpx_action_t deliverSpikes; ///> delivering of (queued) synapses (runs NET_RECEIVE on mod files)

    ///queue of incoming synapses (delivered at the end of the step)
    std::priority_queue<Synapse> queuedSynapses;

  private:

    hpx_t mutex;			///> mutex to protect this branch's memory access

    //For recursive methods
    hpx_t * futures;
    int * futuresSizes;
    void** futuresAddrs;
    FwSubFutureData * futuresData;


    static int setupMatrixRHS_handler(const char isSoma, const double v_parent);
    static int setupMatrixLHS_handler(const char isSoma);
    static int setupMatrixInitValues_handler();
    static int setV_handler(const double v);
    static int callMechsFunction_handler(const Mechanism::Function functionId);
    static int queueSpike_handler(const Synapse * syn, size_t);
    static int deliverSpikes_handler(const int time);
    static int init_handler(const int n, const double *a, const double *b, const double *d,
                            const double *v, const double *rhs, const double *area,
                            const int m, const int * mechsOffsets, const double *data,
                            const Datum *pdata, const int childrenCount, const hpx_t * children);
};

//TODO we're not using this
class FwSubFutureData
{
  public:
    double rhs;
    double d;
} ;

}; //namespace
