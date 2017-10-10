#include "neurox/neurox.h"
#include "stdio.h"
#include "string.h"
#include "tclap/CmdLine.h"

using namespace neurox::tools;

CmdLineParser::CmdLineParser()
    :  // from nrnoptarg.cpp::cn_parameters():
      tstart_(0),
      tstop_(100),
      dt_(0.025),
      dt_io_(0.1),
      celsius_(34),
      voltage_(-65),
      forwardSkip_(0),
      prcellgid_(-1) {
  // TODO: missing some inits
  memset(pattern_stim_, '\0', 512);
  memset(input_path_, '\0', 512);
  memset(output_path_, '\0', 512);
}

CmdLineParser::~CmdLineParser() {}

CmdLineParser::CmdLineParser(int argc, char** argv) : CmdLineParser() {
  Parse(argc, argv);
}

void CmdLineParser::Parse(int argc, char** argv) {
  try {
#ifndef VERSION_STRING
    TCLAP::CmdLine cmd("neurox simulator.", ' ');
#else
    TCLAP::CmdLine cmd("neurox simulator.", ' ', VERSION_STRING);
#endif

    // neurox only command line arguments
    //(NOTE: SwitchArg does not require cmd.add())

    TCLAP::SwitchArg mechs_parallelism(
        "M", "multimex", "activates graph-based parallelism of mechanisms.",
        cmd, false);
    TCLAP::SwitchArg allreduce_at_locality(
        "R", "reduce-by-locality",
        "perform HPX all-reduce operation at "
        "locality level instead of neuron "
        "level (better for small cluster).",
        cmd, false);
    TCLAP::SwitchArg output_statistics(
        "S", "output-statistics",
        "outputs files with memory consumption and mechanism distribution.",
        cmd, false);
    TCLAP::SwitchArg output_mechanisms_dot(
        "1", "output-mechs",
        "outputs mechanisms.dot with mechanisms dependencies.", cmd, false);
    TCLAP::SwitchArg output_netcons_dot(
        "2", "output-netcons",
        "outputs netcons.dot with netcons information across neurons.", cmd,
        false);
    TCLAP::SwitchArg output_compartments_dot(
        "3", "output-compartments",
        "outputs compartments_*.dot files displaying neurons morpholgies.", cmd,
        false);
    TCLAP::SwitchArg load_balancing(
        "L", "load-balancing",
        "performs dynamic load balancing of neurons and branches.", cmd, false);
    TCLAP::ValueArg<int> branch_parallelism_depth(
        "B", "branching-depth",
        "Depth of branches parallelism (0: none, default)", false, 0, "int");
    TCLAP::ValueArg<int> synchronizer(
        "A", "synchronizer",
        "\
[0] BackwardEulerCoreneuronDebug\
\n[1] BackwardEulerWithAllReduceBarrier (default)\
\n[2] BackwardEulerWithSlidingTimeWindow\
\n[3] BackwardEulerWithTimeDependencyLCO\
\n[4] BackwardEulerCoreneuron\
\n[9] All methods sequentially (NOTE: neurons data does not reset)",
        false, (int)synchronizers::Synchronizers::kAllReduce, "int");

    TCLAP::ValueArg<int> interpolator(
        "I", "interpolator",
        "\
[0] CVODES with Diagonal Jacobian solver\
\n[1] CVODES with Dense Jacobian\
\n[2] CVODES with Diagonal Jacobian\
\n[3] CVODES with Sparse Jacobian\
\n[9] Backward Euler (default)",
        false, (int)interpolators::Interpolators::kBackwardEuler, "int");

    cmd.add(branch_parallelism_depth);
    cmd.add(synchronizer);
    cmd.add(interpolator);

    // coreneuron command line parameters
    TCLAP::ValueArg<floble_t> tstart(
        "s", "tstart", "Execution start time (msecs). The default value is 0",
        false, 0, "floble_t");
    cmd.add(tstart);
    TCLAP::ValueArg<floble_t> tstop(
        "e", "tstop", "Execution stop time (msecs). The default value is 10",
        false, 10, "floble_t");
    cmd.add(tstop);
    TCLAP::ValueArg<floble_t> dt(
        "t", "dt",
        "Fixed (or minimum) time step size for fixed (or variable) step interpolation:\
\n - CVODES with Diagonal Jacobian solver: default 0.0001 msecs\
\n - CVODES with Dense Jacobian: default 0.001 msecs\
\n - CVODES with Diagonal Jacobian solver: default 0.00001 msecs\
\n - CVODES with Sparse Jacobian solver: default 0.0001 msecs\
\n - Backward Euler: default 0.025",
        false, DEF_dt, "floble_t");
    cmd.add(dt);
    TCLAP::ValueArg<floble_t> dt_io(
        "i", "dt_io", "I/O time step (msecs). The default value is 0.1", false,
        0.1, "floble_t");
    cmd.add(dt_io);
    TCLAP::ValueArg<floble_t> celsius(
        "l", "celsius",
        "System temperatura (celsius degrees). The default value is 34", false,
        34.0, "floble_t");
    cmd.add(celsius);
    TCLAP::ValueArg<floble_t> forwardskip(
        "k", "forwardskip",
        "Set forwardskip time step (msecs). The default value is 0", false, 0,
        "floble_t");
    cmd.add(forwardskip);
    TCLAP::ValueArg<neuron_id_t> prcellgid(
        "g", "prcellgid",
        "Output prcellstate information for given gid. The default value is -1",
        false, -1, "neuron_id_t");
    cmd.add(prcellgid);
    TCLAP::ValueArg<std::string> pattern_stim(
        "p", "pattern",
        "Apply patternstim with the spike file. No default value", false, "",
        "string");
    cmd.add(pattern_stim);
    TCLAP::ValueArg<std::string> output_path(
        "o", "outputpath",
        "Path to output directory. The default value is ./output", false,
        "./output", "string");
    cmd.add(output_path);
    TCLAP::ValueArg<std::string> input_path("d", "inputpath",
                                            "Path to input files directory",
                                            true, "./input", "string");
    cmd.add(input_path);

    // parse command line arguments
    cmd.parse(argc, argv);

    // copy parameters
    sprintf(this->input_path_, "%s", input_path.getValue().c_str());
    sprintf(this->output_path_, "%s", output_path.getValue().c_str());
    sprintf(this->pattern_stim_, "%s", pattern_stim.getValue().c_str());
    this->tstart_ = tstart.getValue();
    this->tstop_ = tstop.getValue();
    this->dt_ = dt.getValue();
    this->dt_io_ = dt_io.getValue();
    this->celsius_ = celsius.getValue();
    this->prcellgid_ = prcellgid.getValue();
    this->forwardSkip_ = forwardskip.getValue();
    this->voltage_ = DEF_vrest;
    this->second_order_ = (char)DEF_secondorder;
    this->rev_dt_ = 1 / dt.getValue();
    this->celsius_ = DEF_celsius;

    this->output_statistics_ = output_statistics.getValue();
    this->output_mechanisms_dot_ = output_mechanisms_dot.getValue();
    this->output_netcons_dot = output_netcons_dot.getValue();
    this->output_compartments_dot_ = output_compartments_dot.getValue();
    this->mechs_parallelism_ = mechs_parallelism.getValue();
    this->allreduce_at_locality_ = allreduce_at_locality.getValue();
    this->load_balancing_ = load_balancing.getValue();
    this->branch_parallelism_depth_ = branch_parallelism_depth.getValue();
    this->synchronizer_ = (synchronizers::Synchronizers)synchronizer.getValue();
    neurox::synchronizer_ = synchronizers::Synchronizer::New(this->synchronizer_);
    this->interpolator_ = (interpolators::Interpolators)interpolator.getValue();

    if (this->branch_parallelism_depth_ < 0)
      throw TCLAP::ArgException("branch parallism depth should be >= 0",
                                "branch-parallelism-depth");

    if (this->tstop_ <= 0)
      throw TCLAP::ArgException(
          "execution time (ms) should be a positive value", "tstop");
    floble_t remainder_tstop_tcomm =
        fmod(this->tstop_, neurox::min_delay_steps_ * this->dt_);

    if (this->branch_parallelism_depth_ > 0 &&
        this->interpolator_ != interpolators::Interpolators::kBackwardEuler)
      throw TCLAP::ArgException(
          "cant run branch-level parallelism with variable-step methods");

    // handling of dt for variable-step interpolations
    if (this->interpolator_ != Interpolators::kBackwardEuler) {
      if (!dt.isSet())  // if not user-provided
      {
        switch (this->interpolator_) {
          case Interpolators::kCvodePreConditionedDiagSolver:
            this->dt_ = 1e-4;
            break;
          case Interpolators::kCvodeDenseMatrix:
            this->dt_ = 1e-3;
            break;
          case Interpolators::kCvodeDiagonalMatrix:
            this->dt_ = 1e-5;
            break;
          default:
            this->dt_ = 1e-4;
        }
      }
    } else  // ... for fixed-step interpolation
    {
      if (this->dt_ <= 0)
        throw TCLAP::ArgException(
            "time-step size (ms) should be a positive value", "dt");

      if (!(remainder_tstop_tcomm < 0.00001 ||
            remainder_tstop_tcomm >
                neurox::min_delay_steps_ * this->dt_ - 0.00001))
        throw TCLAP::ArgException(
            "execution time " + to_string(this->tstop_) +
                "ms should be a multiple of the communication delay " +
                to_string(neurox::min_delay_steps_ * this->dt_) + " ms",
            "tstop");
    }

  } catch (TCLAP::ArgException& e) {
    printf("TCLAP error: %s (%s).\n", e.error().c_str(), e.argId().c_str());
    exit(1);
  }
}
