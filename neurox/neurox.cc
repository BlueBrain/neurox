#include "neurox/neurox.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <map>

using namespace neurox::synchronizers;
using namespace neurox::interpolators;
using namespace neurox::wrappers;

namespace neurox {

// TODO compute at runtime
double min_synaptic_delay_ = 0.1;
hpx_t *neurons_ =
    nullptr;  // TODO can we get rid of this? replace by hpx_t array
int neurons_count_ = 0;
int mechanisms_count_ = -1;
int *mechanisms_map_ = nullptr;
neurox::Mechanism **mechanisms_ = nullptr;
neurox::tools::CmdLineParser *input_params_ = nullptr;
neurox::synchronizers::Synchronizer *synchronizer_ = nullptr;

Mechanism *GetMechanismFromType(int type) {
  assert(mechanisms_map_[type] != -1);
  return mechanisms_[mechanisms_map_[type]];
}

hpx_action_t Main = 0;
static int Main_handler() {
  printf("\nneurox::Main (localities: %d, threads/locality: %d, %s)\n",
         NumRanks(), NumThreads(), LAYOUT == 0 ? "SoA" : "AoS");
  DebugMessage("neurox::Input::DataLoader::Init...\n");
  CallAllLocalities(neurox::input::DataLoader::Init);
  DebugMessage("neurox::Input::DataLoader::InitMechanisms...\n");
  CallAllLocalities(neurox::input::DataLoader::InitMechanisms);
  DebugMessage("neurox::Input::DataLoader::InitNeurons...\n");
  CallAllLocalities(neurox::input::DataLoader::InitNeurons);
  DebugMessage("neurox::Input::DataLoader::InitNetcons...\n");
  CallAllNeurons(neurox::input::DataLoader::InitNetcons);
  DebugMessage("neurox::Input::DataLoader::Finalize...\n");
  CallAllLocalities(neurox::input::DataLoader::Finalize);
  DebugMessage("neurox::Branch::BranchTree::InitLCOs...\n");
  CallAllNeurons(Branch::BranchTree::InitLCOs);
  neurox::input::Debugger::CompareMechanismsFunctions();
  neurox::input::Debugger::CompareAllBranches();

  if (neurox::input_params_->output_statistics_) {
    tools::Statistics::OutputMechanismsDistribution();
    tools::Statistics::OutputSimulationSize();
    // hpx_exit(0,NULL);
  }

  // call init action on each neuron (e.g. Finitialize, Cvodes init)
  DebugMessage("neurox::Branch::Initialize...\n");
  CallAllNeurons(Branch::Initialize);
#ifndef NDEBUG
  hpx_bcast_rsync(neurox::input::Debugger::Finitialize);
  hpx_bcast_rsync(neurox::input::Debugger::ThreadTableCheck);
  neurox::input::Debugger::CompareAllBranches();
#endif

  // iterator through all synchronizers (if many) and run
  hpx_time_t total_time_now = hpx_time_now();
  const int synchronizer = (int)input_params_->synchronizer_;
  const bool run_all = synchronizer == (int)Synchronizers::kBenchmarkAll;
  const int init_type = run_all ? 0 : synchronizer;
  const int end_type =
      run_all ? (int)Synchronizers::kSynchronizersCount : synchronizer + 1;
  const double tstop = input_params_->tstop_;
  for (int type = init_type; type <= end_type; type++) {
    CallAllLocalities(Synchronizer::InitializeLocality, &type, sizeof(type));
    CallAllNeurons(Synchronizer::InitializeNeuron);

    hpx_time_t time_now = hpx_time_now();
    if (input_params_->locality_comm_reduce_)
      CallAllLocalities(Synchronizer::RunLocality, &tstop, sizeof(tstop));
    else
      CallAllNeurons(Synchronizer::RunNeuron, &tstop, sizeof(tstop));
    double time_elapsed = hpx_time_elapsed_ms(time_now) / 1e3;

    printf("neurox::%s (%d neurons, t=%.03f secs, dt=%.03f milisecs\n",
           synchronizer_->GetString(), neurox::neurons_count_,
           input_params_->tstop_ / 1000, input_params_->dt_);

#ifdef NDEBUG
    // output benchmark info
    printf("csv,%d,%d,%d,%.1f,%.1f,%d,%d,%d,%d,%.2f\n", neurox::neurons_count_,
           hpx_get_num_ranks(), hpx_get_num_threads(),
           neurox::neurons_count_ / (double)hpx_get_num_ranks(),
           input_params_->tstop_, sync_synchronizer_->GetType(),
           input_params_->mechs_parallelism_ ? 1 : 0,
           input_params_->branch_parallelism_depth_,
           input_params_->allreduce_at_locality_ ? 1 : 0, time_elapsed);
    fflush(stdout);
#endif
    CallAllLocalities(Synchronizer::ClearLocality);
    delete synchronizer_;
  }

  double total_elapsed_time = hpx_time_elapsed_ms(total_time_now) / 1e3;
  printf(
      "neurox::end (%d neurons, biological time: %.3f secs, solver time: %.3f "
      "secs).\n",
      neurox::neurons_count_, input_params_->tstop_ / 1000.0,
      total_elapsed_time);

  CallAllNeurons(Branch::Clear);
  hpx_bcast_rsync(neurox::Clear);
  hpx_exit(0, NULL);
}

hpx_action_t Clear = 0;
int Clear_handler() {
  NEUROX_MEM_PIN(uint64_t);
  delete[] neurox::mechanisms_;
  delete[] neurox::neurons_;
  delete[] neurox::mechanisms_map_;
  delete synchronizer_;

#ifndef NDEBUG
  input::DataLoader::CleanCoreneuronData(true);
#endif
  return MemoryUnpin(target);
}

void DebugMessage(const char *str) {
#ifndef NDEBUG
  printf("%s", str);
  fflush(stdout);
#endif
}

bool ParallelExecution() { return hpx_get_num_ranks() > 1; }

void RegisterHpxActions() {
  RegisterZeroVarAction(Main, Main_handler);
  RegisterZeroVarAction(Clear, Clear_handler);
}

};  // neurox
