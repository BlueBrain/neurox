#include "neurox/neurox.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <map>

using namespace neurox::synchronizers;
using namespace neurox::interpolators;

namespace neurox {

//TODO compute at runtime
double min_synaptic_delay_ = 0.1 + 0.00000001;
hpx_t *neurons_ = nullptr;
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
         neurox::wrappers::NumRanks(), neurox::wrappers::NumThreads(),
         LAYOUT == 0 ? "SoA" : "AoS");
  DebugMessage("neurox::Input::DataLoader::Init...\n");
  neurox::wrappers::CallAllLocalities(neurox::input::DataLoader::Init);
  DebugMessage("neurox::Input::DataLoader::InitMechanisms...\n");
  neurox::wrappers::CallAllLocalities(
      neurox::input::DataLoader::InitMechanisms);
  DebugMessage("neurox::Input::DataLoader::InitNeurons...\n");
  neurox::wrappers::CallAllLocalities(neurox::input::DataLoader::InitNeurons);
  DebugMessage("neurox::Input::DataLoader::InitNetcons...\n");
  neurox::wrappers::CallAllNeurons(neurox::input::DataLoader::InitNetcons);
  DebugMessage("neurox::Input::DataLoader::Finalize...\n");
  neurox::wrappers::CallAllLocalities(neurox::input::DataLoader::Finalize);
  DebugMessage("neurox::Branch::BranchTree::InitLCOs...\n");
  neurox::wrappers::CallAllNeurons(Branch::BranchTree::InitLCOs);
  neurox::input::Debugger::CompareMechanismsFunctions();
  neurox::input::Debugger::CompareAllBranches();

  if (neurox::input_params_->output_statistics_) {
    tools::Statistics::OutputMechanismsDistribution();
    tools::Statistics::OutputSimulationSize();
    // hpx_exit(0,NULL);
  }

  //call init action on each neuron (e.g. Finitialize, Cvodes init)
  DebugMessage("neurox::Branch::Initialize...\n");
  neurox::wrappers::CallAllNeurons(Branch::Initialize);
#ifndef NDEBUG
  hpx_bcast_rsync(neurox::input::Debugger::Finitialize);
  hpx_bcast_rsync(neurox::input::Debugger::ThreadTableCheck);
  neurox::input::Debugger::CompareAllBranches();
#endif

  //iterator through all synchronizers (if many) and run
  hpx_time_t total_time_now = hpx_time_now();
  const int synchronizer = (int) input_params_->synchronizer_;
  const bool run_all = synchronizer == (int) Synchronizers::kBenchmarkAll;
  const int init_type = run_all ? 0 : synchronizer;
  const int end_type = run_all ? (int) Synchronizers::kSynchronizersCount : synchronizer;
  const double tstop = input_params_->tstop_;
  for (int type = init_type; type < end_type; type++)
  {
      wrappers::CallAllLocalities(Synchronizer::InitLocality, &type, sizeof(type));

      hpx_time_t time_now = hpx_time_now();
      if (input_params_->locality_comm_reduce_)
        neurox::wrappers::CallAllLocalities(Synchronizer::RunLocality, &tstop, sizeof(tstop));
      else
        neurox::wrappers::CallAllNeurons(Synchronizer::RunNeuron, &tstop, sizeof(tstop));
      double time_elapsed = hpx_time_elapsed_ms(time_now) / 1e3;

      printf("neurox::%s (%d neurons, t=%.03f secs, dt=%.03f milisecs\n",
            synchronizer_->GetString(), neurox::neurons_count_, input_params_->tstop_ / 1000,
            input_params_->dt_);

#ifdef NDEBUG
      // output benchmark info
      printf("csv,%d,%d,%d,%.1f,%.1f,%d,%d,%d,%d,%.2f\n",
             neurox::neurons_count_, hpx_get_num_ranks(), hpx_get_num_threads(),
             neurox::neurons_count_ / (double)hpx_get_num_ranks(),
             input_params_->tstop_, sync_synchronizer_->GetType(),
             input_params_->mechs_parallelism_ ? 1 : 0,
             input_params_->branch_parallelism_depth_,
             input_params_->allreduce_at_locality_ ? 1 : 0, time_elapsed);
      fflush(stdout);
#endif
      neurox::wrappers::CallAllLocalities(Synchronizer::ClearLocality);
      delete synchronizer_;
  }

  DebugMessage("neurox::Interpolator::ClearNeuron...\n");
  wrappers::CallAllNeurons(Branch::Clear);
  hpx_bcast_rsync(neurox::Clear);

  double total_elapsed_time = hpx_time_elapsed_ms(total_time_now) / 1e3;
  printf("neurox::end (%d neurons, biological time: %.3f secs, solver time: %.3f "
      "secs).\n",
      neurox::neurons_count_, input_params_->tstop_ / 1000.0, total_elapsed_time);
  hpx_exit(0, NULL);
}

hpx_action_t Clear = 0;
int Clear_handler() {
  NEUROX_MEM_PIN(uint64_t);
  delete[] neurox::mechanisms_;
  delete[] neurox::neurons_;
  delete[] neurox::mechanisms_map_;

  if (input_params_->locality_comm_reduce_) {
    AllreduceSynchronizer::AllReducesInfo::AllReduceLocality::locality_neurons_
        ->clear();
    delete AllreduceSynchronizer::AllReducesInfo::AllReduceLocality::
        locality_neurons_;
    AllreduceSynchronizer::AllReducesInfo::AllReduceLocality::
        locality_neurons_ = nullptr;
  }

#ifndef NDEBUG
  input::DataLoader::CleanCoreneuronData(true);
#endif
  return wrappers::MemoryUnpin(target);
}

void DebugMessage(const char *str) {
#ifndef NDEBUG
  printf("%s", str);
  fflush(stdout);
#endif
}

bool ParallelExecution() { return hpx_get_num_ranks() > 1; }

void RegisterHpxActions() {
  wrappers::RegisterZeroVarAction(Main, Main_handler);
  wrappers::RegisterZeroVarAction(Clear, Clear_handler);
}

};  // neurox
