#pragma once

#include "hpx/hpx.h"

#include "coreneuron/nrnoc/membfunc.h"

typedef hpx_addr_t hpx_t;
typedef unsigned char byte;

class GlobalInfo;
class Neuron;
class Branch;
class FwSubFutureData;

/**
 * @brief The GlobalInfo class
 * Represents all information in the system including execution parameters, neurons and synapses
 */
class GlobalInfo
{
  public:

    GlobalInfo();
    ~GlobalInfo();

    //global vars: all localities hold the same value
    int neuronsCount; 	///> total neurons count in the system
    hpx_t neuronsAddr; 	///> hpx address of the first position of the neurons array
    int multiSplit; 	///> 0 or 1 for multisplit or not

    //Execution parameters (cn_input_parameters)
    int secondorder; 	///> 0 means crank-nicolson. 2 means currents adjusted to t+dt/2
    double t; 			///> current simulation time (msecs)
    double dt; 			///> time step ie delta-t (msecs)
    double rev_dt; 		///> reverse of delta t (1/msecs)
    double celsius; 	///> celsius temperature (degrees)
    double tstart; 		///> start time of simulation in msec*/
    double tstop;		///> stop time of simulation in msec*/
    double dt_io;    	///> i/o timestep to use in msec*/
    double voltage;     ///> TODO: what's this?
    double maxdelay;    ///> TODO: do we need this?
    double mindelay;    ///> minimum synaptic delay
    double forwardSkip;	///> forward skip time
    int prcellgid;      ///> gid of cell for prcellstate

    char inputPath[2048];		///>path of input directory
    char outputPath[2048];		///>path of output directory
    char patternStimFile[2048];	///>patternStim file path
} ;

/**
 * @brief The Neuron class
 * Represents a neuron as its metadata and a tree-based orphology
 */
class Neuron
{
  public:
    hpx_t topBranch;		///> hpx address of the top compartment (soma)

    //neuron metadata
    int id;					///> neuron global id
};

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
};

class FwSubFutureData
{
  public:
    double rhs;
    double d;
} ;
