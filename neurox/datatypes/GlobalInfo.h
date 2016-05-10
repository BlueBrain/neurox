#pragma once

#include "neurox/neurox.h"

/**
 * @brief The GlobalInfo class
 * Represents all information in the system including execution parameters, neurons and synapses
 */
class GlobalInfo
{
  public:

    // default values from register_mech.c:: initnrn()
    GlobalInfo();
    GlobalInfo(int argc, char ** argv);
    ~GlobalInfo();

    //global vars: all localities hold the same value
    int neuronsCount; 	///> total neurons count in the system
    int multiSplit; 	///> 0 or 1 for multisplit or not
    hpx_t neuronsAddr; 	///> hpx address of the first position of the neurons array

    //Execution parameters (cn_input_parameters)
    int secondorder; 	///> 0 means crank-nicolson. 2 means currents adjusted to t+dt/2
    int prcellgid;      ///> gid of cell for prcellstate
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

    char inputPath[2048];	///>path of input directory
    char outputPath[2048];	///>path of output directory
    char patternStim[2048];	///>patternStim file path

  private:
    void parseCommandLine(int argc, char ** argv); /// populates values from cmd line
} ;
