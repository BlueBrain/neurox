#pragma once

#include "neurox/neurox.h"

namespace neurox
{

namespace Input
{
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
    double t; 			///> current simulation time (msecs)
    double dt; 			///> time step ie delta-t (msecs)
    double rev_dt; 		///> reverse of delta t (1/msecs)
    double celsius; 	///> celsius temperature (degrees)
    double tstart; 		///> start time of simulation in msec*/
    double tstop;		///> stop time of simulation in msec*/
    double dt_io;    	///> i/o timestep to use in msec*/
    double voltage;     ///> initial voltage set on all neurons
    double forwardSkip;	///> forward skip time

    char inputPath[512];	///> path of input directory
    char outputPath[512];	///> path of output directory
    char patternStim[512];	///> patternStim file path (the filename of an output_spikes.h format raster file.)

  private:
    /// Parses command line arguments and populates structure
    void parseCommandLine(int argc, char ** argv);

};
}; //InputParams
}; //NeuroX

