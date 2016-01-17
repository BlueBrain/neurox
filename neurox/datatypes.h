#pragma once

#include "hpx/hpx.h"

typedef hpx_addr_t hpx_t;

typedef struct fwSubFutureData
{
    double rhs;
    double d;
} FwSubFutureData;

class Neuron
{
  public:
    hpx_t branches;
    int branchesCount;
    int neuronMetaData;
};

class Branch
{
  public:
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
    
  private:
    //semaphore
    hpx_t mutex;
    
#if USE_LCO_FUTURE_ARRAY==0
    //For fwSub Method
    hpx_t * futures;
    int * futuresSizes;
    void** futuresAddrs;
    FwSubFutureData * futuresData;
#endif    
};

class GlobalInfo
{
  public:

    GlobalInfo():
        branchesCount(0), compartmentsCount(0), neuronsCount(0){}

    //global vars: all localities hold the same value
    int neuronsCount;
    hpx_t neuronsAddr;
    long long compartmentsCount;
    long long branchesCount;
    int multiSplit;

    //This mutex forces only one hdf5 to be opened at a time. localities hold a diff value
    //hpx_t mutex_io_h5;

} ;
