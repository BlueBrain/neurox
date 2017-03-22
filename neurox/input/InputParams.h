#pragma once
#include "neurox.h"

namespace neurox
{

enum Algorithm {
    BackwardEulerSyncFixedCommStepDebug=0,
    BackwardEulerAsyncFixedCommStep=1,
    BackwardEulerWithPairwiseSteping=2
};

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
    char secondorder;       ///> 0 means crank-nicolson. 2 means currents adjusted to t+dt/2
    neuron_id_t prcellgid;  ///> gid of cell for prcellstate
    floble_t dt; 	        ///> time step ie delta-t (msecs)
    floble_t rev_dt; 	    ///> reverse of delta t (1/msecs)
    floble_t celsius; 	    ///> celsius temperature (degrees)
    floble_t tstart; 	    ///> start time of simulation in msec*/
    floble_t tstop;	        ///> stop time of simulation in msec*/
    floble_t dt_io;    	    ///> i/o timestep to use in msec*/
    floble_t voltage;       ///> initial voltage set on all neurons
    floble_t forwardSkip;   ///> forward skip time

    char inputPath[512];    ///> path of input directory
    char outputPath[512];   ///> path of output directory
    char patternStim[512];  ///> patternStim file path (the filename of an output_spikes.h format raster file.)

    //neurox specific options
    bool outputStatistics;      ///> outputs statistics file
    bool outputMechanismsDot;   ///> outputs mechanisms.dot file
    bool outputNetconsDot;      ///> outputs netcons.dot file
    bool outputCompartmentsDot; ///> outputs compartments*.dot files
    bool multiMex;              ///> graph-based parallelism of mechanisms
    bool multiSplix;            ///> tree-based parallelism of morphologies

    bool parallelDataLoading;///>coreneuron mpi-based execution

    neurox::Algorithm algorithm;

  private:
    /// Parses command line arguments and populates structure
    void parseCommandLine(int argc, char ** argv);
};
}; //InputParams
}; //NeuroX

