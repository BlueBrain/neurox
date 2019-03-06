#include "neurox/neurox.h"
#include "stdio.h"
#include "string.h"
#include "tclap/CmdLine.h"

using namespace neurox::tools;
using namespace neurox::interpolators;
using namespace neurox::synchronizers;

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
    TCLAP::SwitchArg linearize_containers(
        "W", "linearize-containers", "use linear represnetation containers",
        cmd, false);
    TCLAP::SwitchArg output_statistics(
        "5", "output-statistics",
        "outputs files with memory consumption and mechanism distribution.",
        cmd, false);
    TCLAP::SwitchArg output_comm_count(
        "4", "output-comm-count",
        "outputs final information about number of point-to-point and reduce "
        "operation.",
        cmd, false);
    TCLAP::SwitchArg output_compartments_dot(
        "3", "output-compartments",
        "outputs compartments_*.dot files displaying neurons morpholgies.", cmd,
        false);
    TCLAP::SwitchArg output_netcons_dot(
        "2", "output-netcons",
        "outputs netcons.dot with netcons information across neurons.", cmd,
        false);
    TCLAP::SwitchArg output_mechanisms_dot(
        "1", "output-mechs",
        "outputs mechanisms.dot with mechanisms dependencies.", cmd, false);

    TCLAP::SwitchArg locality_comm_reduce("C", "comm-reduce",
                                          "perform HPX all-reduce operation at "
                                          "locality level instead of neuron "
                                          "level (better for small cluster).",
                                          cmd, false);

    TCLAP::SwitchArg neurons_scheduler(
        "S", "scheduler", "last neuron goes first scheduler", cmd, false);

    TCLAP::SwitchArg graph_mechs_parallelism(
        "G", "graph-parallelism",
        "activates graph-based parallelism of mechanisms.", cmd, false);

    TCLAP::ValueArg<floble_t> subsection_complexity(
        "", "L",
        "scale constant to groups of subsection complexity (constant k' for "
        "branch parallelism on multiple localities)",
        false, 0, "floble_t");
    cmd.add(subsection_complexity);

    TCLAP::SwitchArg load_balancing(
        "L", "load-balancing",
        "performs dynamic load balancing of neurons and branches.", cmd, false);

    TCLAP::ValueArg<floble_t> subtree_complexity(
        "", "B",
        "scale constant to subsection complexity (constant k for branch "
        "parallelism)",
        false, 0.3, "floble_t");
    cmd.add(subtree_complexity);

    TCLAP::SwitchArg branch_parallelism(
        "B", "branch-parallelism",
        "performs branch-level parallelism on neurons", cmd, false);

    TCLAP::ValueArg<floble_t> mech_instance_percent_per_block(
        "", "M",
        "percentage of total workload assigned to each mechanism instances "
        "block (for mech-instance parallelism)",
        false, 0.1, "floble_t");
    cmd.add(mech_instance_percent_per_block);

    TCLAP::ValueArg<int> processor_cache_line_size_l1(
        "", "L1", "processor Level-1 cache line size", false, 64, "int");
    cmd.add(processor_cache_line_size_l1);

    TCLAP::SwitchArg mech_instances_parallelism(
        "M", "mech-parallelism", "parallelism of mechanisms instances", cmd,
        false);

    TCLAP::ValueArg<int> synchronizer(
        "A", "synchronizer",
        "\
\n[0] Time Dependency LCO\
\n[1] All-reduce barrier (default)\
\n[2] Sliding Time Window\
\n[3] MPI-based (a la Coreneuron)\
\n[8] Sequential Single-step Barrier (debug  only)",
        false, (int)synchronizers::SynchronizerIds::kAllReduce, "int");
    cmd.add(synchronizer);

    TCLAP::ValueArg<int> interpolator(
        "I", "interpolator",
        "\
[0] CVODE with Preconditioned Diagonal Jacobian (a la NEURON)\
\n[1] CVODE with Dense Jacobian\
\n[2] CVODE with Diagonal Jacobian\
\n[3] CVODE with Sparse Jacobian\
\n[9] Backward Euler (default)",
        false, (int)interpolators::InterpolatorIds::kBackwardEuler, "int");
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
\n - Varible timestepping default (NEURON value) is 0:\
\n   - CVODE with Preconditioned Diagonal Jacobian solver: recommended 0.0001 msecs\
\n   - CVODE with Dense Jacobian: recommended 0.001 msecs\
\n   - CVODE with Diagonal Jacobian solver: recommended 0.00001 msecs\
\n   - CVODE with Sparse Jacobian solver: recommended 0.0001 msecs\
\n - Backward Euler: default (NEURON value) 0.025",
        false, DEF_dt, "floble_t");
    cmd.add(dt);

    /* From CVODE 4.0.2 manual:
     * The scalar relative tolerance reltol is to be set to control relative
     * errors. So reltol=1e−4 means that errors are controlled to .01%. We do
     * not recommend using reltol larger than 1e−3. On the other hand, reltol
     * should not be so small that it is comparable to the unit roundoff of the
     * machine arithmetic (generally around 1.0E-15). */
    TCLAP::ValueArg<floble_t> cvode_rtol(
        "", "rtol",
        "relative tolerance for variable timestepping. Default value is 0 "
        "(a la NEURON, not recommended by CVODE manual). Recommended: 1e-4.",
        false, 0, "floble_t");
    cmd.add(cvode_rtol);

    TCLAP::ValueArg<floble_t> cvode_atol(
        "", "atol",
        "absolute tolerance for voltages and states in variable timestepping. "
        "Default (NEURON) value is 1e-3.",
        false, 1e-3, "floble_t");
    cmd.add(cvode_atol);

    TCLAP::ValueArg<floble_t> cvode_event_group(
        "", "queue_group",
        "interval (msecs) for grouping of events in variable timestepping. "
        "Default (NEURON) value is 0.",
        false, 0, "floble_t");
    cmd.add(cvode_event_group);

    TCLAP::SwitchArg cvode_speculative(
        "", "cvode-speculate", "perform CVODE speculative stepping", cmd,
        false);

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

    this->linearize_containers_ = linearize_containers.getValue();
    this->output_statistics_ = output_statistics.getValue();
    this->output_comm_count_ = output_comm_count.getValue();
    this->output_mechanisms_dot_ = output_mechanisms_dot.getValue();
    this->output_netcons_dot = output_netcons_dot.getValue();
    this->output_compartments_dot_ = output_compartments_dot.getValue();
    this->graph_mechs_parallelism_ = graph_mechs_parallelism.getValue();
    this->mech_instances_parallelism_ = mech_instances_parallelism.getValue();
    this->locality_comm_reduce_ = locality_comm_reduce.getValue();
    this->scheduler_ = neurons_scheduler.getValue();
    this->load_balancing_ = load_balancing.getValue();
    this->branch_parallelism_ = branch_parallelism.getValue();
    this->synchronizer_ =
        (synchronizers::SynchronizerIds)synchronizer.getValue();
    neurox::synchronizer_ =
        synchronizers::Synchronizer::New(this->synchronizer_);
    this->interpolator_ =
        (interpolators::InterpolatorIds)interpolator.getValue();

    this->subtree_complexity = subtree_complexity.getValue();
    this->subsection_complexity = subsection_complexity.getValue();
    this->mech_instance_percent_per_block =
        mech_instance_percent_per_block.getValue();
    this->processor_cache_line_size_l1_ =
        processor_cache_line_size_l1.getValue();

    if (this->tstop_ <= 0)
      throw TCLAP::ArgException(
          "execution time (ms) should be a positive value", "tstop");
    floble_t remainder_tstop_tcomm =
        fmod(this->tstop_, neurox::min_synaptic_delay_);

    if (this->branch_parallelism_ &&
        this->interpolator_ != interpolators::InterpolatorIds::kBackwardEuler)
      throw TCLAP::ArgException(
          "cant run branch-level parallelism with variable-step methods");

    if (!(fabs(remainder_tstop_tcomm) > 0.000000001 ||
          fabs(remainder_tstop_tcomm) <
              neurox::min_synaptic_delay_ - 0.000000001))
      throw TCLAP::ArgException(
          "execution time " + to_string(this->tstop_) +
              "ms should be a multiple of the communication delay " +
              to_string(neurox::min_synaptic_delay_) + " ms",
          "tstop");

    // handling of default dt for variable-step interpolations
    if (this->interpolator_ != InterpolatorIds::kBackwardEuler) {
      this->cvode_rtol_ = cvode_rtol.getValue();
      this->cvode_atol_v_ = cvode_atol.getValue();
      this->cvode_atol_states_ = cvode_atol.getValue();
      this->cvode_event_group_ = cvode_event_group.getValue();
      this->cvode_speculative_ = cvode_speculative.getValue();
      if (!dt.isSet())  // if not user-provided
        this->dt_ = 0;
      if (this->dt_ < 0)
        throw TCLAP::ArgException(
            "time-step size (ms) should be a positive value", "dt");
    }
  } catch (TCLAP::ArgException& e) {
    printf("TCLAP error: %s (%s).\n", e.error().c_str(), e.argId().c_str());
    exit(1);
  }
}
