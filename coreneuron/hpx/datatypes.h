#pragma once

#include "hpx/hpx.h"

typedef hpx_addr_t hpx_t;

void convert_from_coreneuron_to_hpx_datatypes();

typedef struct fwSubFutureData
{
    double rhs;
    double d;
} FwSubFutureData;

typedef struct neuron
{
    hpx_t topBranch;
    int neuronMetaData;
} Neuron;

typedef struct branch
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
    hpx_t * children;
    
    //semaphore
    hpx_t mutex;
    
#if USE_LCO_FUTURE_ARRAY==0
    //For fwSub Method
    hpx_t * futures;
    int * futuresSizes;
    void** futuresAddrs;
    FwSubFutureData * futuresData;
#endif    
} Branch;

typedef struct globalVarsStruct
{
    //global vars: all localities hold the same value
    int neuronsCount;
    hpx_t neuronsAddr;
    int totalCompartmentsCount;
    int totalBranchesCount;

    //This mutex forces only one hdf5 to be opened at a time. localities hold a diff value
    hpx_t mutex_io_h5;

    struct settings
    {
      int multiSplit;
      int neuronsCount;
    } Settings;

} GlobalVars;

