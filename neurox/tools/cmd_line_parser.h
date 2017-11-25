#pragma once
#include "neurox/neurox.h"

namespace neurox {

// Forward declarations
namespace synchronizers {
enum class SynchronizerIds : int;
class Synchronizer;
}

namespace interpolators {
enum class InterpolatorIds : int;
}

namespace tools {
/**
 * @brief The CmdLineParser class
 * Represents all input parameters
 */
class CmdLineParser {
 public:
  // default values from register_mech.c:: initnrn()
  CmdLineParser();
  CmdLineParser(int argc, char** argv);
  ~CmdLineParser();

  // Execution parameters (cn_input_parameters)
  char second_order_;  ///> 0 means crank-nicolson. 2 means currents adjusted to
                       /// t+dt/2
  neuron_id_t prcellgid_;  ///> gid of cell for prcellstate
  floble_t dt_;            ///> time step ie delta-t (msecs)
  floble_t rev_dt_;        ///> reverse of delta t (1/msecs)
  floble_t celsius_;       ///> celsius temperature (degrees)
  floble_t tstart_;        ///> start time of simulation in msec*/
  floble_t tstop_;         ///> stop time of simulation in msec*/
  floble_t dt_io_;         ///> i/o timestep to use in msec*/
  floble_t voltage_;       ///> initial voltage set on all neurons
  floble_t forwardSkip_;   ///> forward skip time

  char input_path_[512];    ///> path of input directory
  char output_path_[512];   ///> path of output directory
  char pattern_stim_[512];  ///> patternStim file path (the filename of an
                            /// output_spikes.h format raster file.)

  // TODO make const?
  // neurox specific options
  bool output_statistics_;        ///> outputs statistics file
  bool output_mechanisms_dot_;    ///> outputs mechanisms.dot file
  bool output_netcons_dot;        ///> outputs netcons.dot file
  bool output_compartments_dot_;  ///> outputs compartments*.dot files
  bool mechs_parallelism_;        ///> graph-based parallelism of mechanisms
  bool locality_comm_reduce_;     ///> locality-based communication reduction

  /// Whether to perform dynamic load balancing of nodes and branches
  bool load_balancing_;

  /// depth tree-based parallelism of morphologies
  bool branch_parallelism_;

  /// neurons sychronization synchronizer
  synchronizers::SynchronizerIds synchronizer_;

  /// interpolation agorithm
  interpolators::InterpolatorIds interpolator_;

 private:
  /// Parses command line arguments and populates structure
  void Parse(int argc, char** argv);
};
};  // CmdLineParser
};  // neurox
