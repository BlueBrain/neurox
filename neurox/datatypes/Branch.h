#pragma once

#include "neurox/neurox.h"
#include "coreneuron/nrnoc/membfunc.h"
#include "coreneuron/nrnoc/membdef.h"

class Branch;
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
    Branch();
    ~Branch();

    void serialize(byte *& bytes_out, int & size_out);
    void deserialize(const byte * bytes_in, const int size_in);

    //nodes
    short int n;			///> number of compartments
    double * b;				///> bottom diagonal of Linear Algebra sparse tridiagonal matrix
    double * d;				///> main diagonal of Linear Algebra spart tridiagonal matrix
    double * a;				///> top diagonal of Linear Algebra sparse tridiagonal matrix
    double * rhs;			///> right-hand side (solution vector) of Linear Algebra solver
    double * v;				///> current voltage per compartment
    double * area;			///> current area per compartment

    //Mechanisms data
    int* is_art;
    int* layout;
    int* is_ion;
    char* pnttype;
    Memb_func** membfunc;
    int* nodesIds;
    int* instanceCount;
    int* dataSize;			///> sizes of double data used per mechanisms
    int* pdataSize;			///> sizes of int data used per mechanisms
    double* data;			///> all double data used by all mechanisms
    int* pdata;				///> all int data used by all mechanisms

    //List of children
    int childrenCount;		///> number of children branches (always >1)
    hpx_t * children;		///> hpx address of the children branches

    static void registerHpxActions(); ///> Register all HPX actions
    static hpx_action_t initialize; ///> Initializes Branch

  private:
    //semaphore
    hpx_t mutex;			///> mutex to protect this branch's memory acces

#if USE_LCO_FUTURE_ARRAY==0
    //For fwSub Method
    hpx_t * futures;
    int * futuresSizes;
    void** futuresAddrs;
    FwSubFutureData * futuresData;
#endif    

    static int initialize_handler(const byte * branch_serial_input,  const size_t size);
};

class FwSubFutureData
{
  public:
    double rhs;
    double d;
} ;
