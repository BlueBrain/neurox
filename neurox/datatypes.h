#pragma once

#include "hpx/hpx.h"

typedef hpx_addr_t hpx_t;
typedef unsigned char byte;

class GlobalInfo;
class Neuron;
class Branch;
class Mechanism;

typedef struct fwSubFutureData
{
    double rhs;
    double d;
} FwSubFutureData;

class Neuron
{
  public:
    hpx_t topBranch;

    //neuron metadata
    int id;
};

class Branch
{
  public:
    Branch();
    ~Branch();

    void serialize(byte *& bytes_out, int & size_out);
    void deserialize(byte * bytes_in, int & size_in);

    //nodes
    short int n;
    double * b;
    double * d;
    double * a;
    double * rhs;
    double * v;
    double * area;

    //Mechanisms data
    int* is_art, layout, is_ion;
    char* pnttype;
    Memb_func** membfunc;
    int* nodesIds, instanceCount, dataSize, pdataSize;
    double* data;
    Datum* pdata;

    //List of children
    int childrenCount;
    hpx_t * children;
    
  private:
    //semaphore
    hpx_t mutex;

    int * pdata;
    double * data;

#if USE_LCO_FUTURE_ARRAY==0
    //For fwSub Method
    hpx_t * futures;
    int * futuresSizes;
    void** futuresAddrs;
    FwSubFutureData * futuresData;
#endif    
};

class Mechanism
{
    int type;
    int * pdata;
    double * data;
    int pdataSize;
    int dataSize;
};

class GlobalInfo
{
  public:

    GlobalInfo():
        neuronsCount(0), neuronsAddr(HPX_NULL), multiSplit(0){}

    //global vars: all localities hold the same value
    int neuronsCount; 	///> total neurons count in the system
    hpx_t neuronsAddr; 	///> hpx address of the first position of the neurons array
    int multiSplit; 	///> 0 or 1 for multisplit or not

    //Execution parameters (cn_input_parameters)
    int secondorder; 	///> 0 means crank-nicolson. 2 means currents adjusted to t+dt/2
    double t; 			///> current simulation time (msecs)
    double dt; 			///> delta t i.e time step (msecs)
    double rev_dt; 		///> reverse of delta t (1/msecs)
    double celsius; 	///> celsius temperature (degrees)
    double tstart; 		///> start time of simulation in msec*/
    double tstop;		///> stop time of simulation in msec*/
    double dt;			///> timestep to use in msec*/
    double dt_io;    	///> i/o timestep to use in msec*/
    double voltage;     ///> TODO: what's this?
    double maxdelay;    ///> TODO: do we need this?
} ;
