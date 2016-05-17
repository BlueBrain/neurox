#pragma once

#include "neurox/neurox.h"

/**
 * @brief The InputParams class
 * Represents all input parameters
 */
class InputParams
{
  public:

    // default values from register_mech.c:: initnrn()
    InputParams();
    InputParams(int argc, char ** argv);
    ~InputParams();

    //Execution parameters (cn_input_parameters)
    int secondorder; 	///> 0 means crank-nicolson. 2 means currents adjusted to t+dt/2
    int prcellgid;      ///> gid of cell for prcellstate
    int multiSplit; 	///> 0 or 1 for multisplit or not
    double t; 			///> current simulation time (msecs)
    double dt; 			///> time step ie delta-t (msecs)
    double rev_dt; 		///> reverse of delta t (1/msecs)
    double celsius; 	///> celsius temperature (degrees)
    double tstart; 		///> start time of simulation in msec*/
    double tstop;		///> stop time of simulation in msec*/
    double dt_io;    	///> i/o timestep to use in msec*/
    double voltage;     ///> initial voltage set on all neurons
    double maxdelay;    ///> TODO: do we need this?
    double mindelay;    ///> minimum synaptic delay
    double forwardSkip;	///> forward skip time

    char inputPath[2048];	///> path of input directory
    char outputPath[2048];	///> path of output directory
    char patternStim[2048];	///> patternStim file path

    static void registerHpxActions();		///> Register all HPX actions
    static hpx_action_t initialize;	///> Initializes InputParams

  private:
    /// Parses command line arguments and populates structure
    void parseCommandLine(int argc, char ** argv);
    static int initialize_handler(const InputParams * inputParams, const size_t size);
} ;

