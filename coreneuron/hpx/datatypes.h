#pragma once

#include "hpx/hpx.h"

typedef struct neuronData
{
    int n;
    double *b;
    double *d;
    double *a;
    double *rhs;
    int *p;
    char path[2048];
    int branchesCount;
    int terminalBranchesCount;
} Neuron_serial;

typedef struct fwSubFutureData
{
    double rhs;
    double d;
} FwSubFutureData;

typedef struct neuron
{
    int solversCount;
    hpx_addr_t solvers_addr;
    int compartmentsCount;
    char path[2048];
} Neuron_GAS;

typedef struct solverData
{   
    short int n;
    short int id;
    double * b;
    double * d;
    double * a;
    double * rhs;
    double * p;
    
    //List of children
    int childrenCount;
    hpx_addr_t * children;
    
    //semaphore
    hpx_addr_t mutex;
    
#if USE_LCO_FUTURE_ARRAY==0
    //For fwSub Method
    hpx_addr_t * futures;	 
    int * sizes;	 
    void** addrs;	 
    FwSubFutureData * futuresData;
#endif    
} Solver;

typedef struct globalVarsStruct
{
    //This mutex forces only one hdf5 to be opened at a time. localities hold a diff value
    hpx_addr_t mutex_io_h5;

    //global vars: all localities hold the same value
    int multiSplit;
    int neuronsCount;
    hpx_addr_t neurons_addr;
    char inputFolder[2048];
} GlobalVars;

