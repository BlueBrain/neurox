/*
# =============================================================================
# Copyright (c) 2015 - 2021-2022 Blue Brain Project/EPFL
#
# See top-level LICENSE file for details.
# =============================================================================.
*/
#include "neurox/neurox.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <map>

using namespace neurox::synchronizers;
using namespace neurox::interpolators;
using namespace neurox::wrappers;

namespace neurox {

hpx_t *neurons_ = nullptr;
int neurons_count_ = 0;
int mechanisms_count_ = -1;
int *mechanisms_map_ = nullptr;
neurox::Mechanism **mechanisms_ = nullptr;
neurox::tools::CmdLineParser *input_params_ = nullptr;
neurox::synchronizers::Synchronizer *synchronizer_ = nullptr;

// locality info
std::vector<hpx_t> *locality::neurons_ = nullptr;
map<neuron_id_t, vector<hpx_t>> *locality::netcons_branches_ = nullptr;
map<neuron_id_t, vector<hpx_t>> *locality::netcons_somas_ = nullptr;
set<pair<floble_t, hpx_t>> *locality::scheduler_neurons_ = nullptr;
libhpx_mutex_t locality::scheduler_lock_ = HPX_NULL;
libhpx_cond_t locality::scheduler_wait_condition_ = HPX_NULL;
hpx_t locality::scheduler_neurons_sema_ = HPX_NULL;
unsigned locality::scheduler_remaining_neurons_ = -1;

#if defined(PRINT_TIME_DEPENDENCY) or defined(PRINT_TIME_DEPENDENCY_MUTEX) or \
    defined(PRINT_TIME_DEPENDENCY_STEP_SIZE)
std::map<hpx_t, neuron_id_t> *locality::from_hpx_to_gid = nullptr;
int locality::scheduler_sema_counter_ = 0;
#endif

Mechanism *GetMechanismFromType(int type) {
  assert(mechanisms_map_[type] != -1);
  return mechanisms_[mechanisms_map_[type]];
}

hpx_action_t Main = 0;
static int Main_handler() {
  hpx_time_t total_time_now = hpx_time_now();
  printf("\nneurox::Main (localities: %d, threads/locality: %d, %s)\n",
         NumRanks(), NumThreads(), LAYOUT == 0 ? "SoA" : "AoS");
  DebugMessage("neurox::input::DataLoader::Init...\n");
  CallAllLocalities(input::DataLoader::Init);
  DebugMessage("neurox::input::DataLoader::InitMechanisms...\n");
  CallAllLocalities(input::DataLoader::InitMechanisms);
  DebugMessage("neurox::input::DataLoader::InitNeurons...\n");
  CallAllLocalities(input::DataLoader::InitNeurons);
  if (input_params_->graph_mechs_parallelism_) {
    DebugMessage("neurox::Branch::InitMechanismsGraph...\n");
    CallAllNeurons(Branch::InitMechanismsGraph);
  }
  if (input_params_->mech_instances_parallelism_) {
    DebugMessage("neurox::Branch::InitMechParallelism...\n");
    CallAllNeurons(Branch::InitMechParallelism);
  }
  DebugMessage("neurox::input::DataLoader::InitNetcons...\n");
  CallAllNeurons(input::DataLoader::InitNetcons);
  DebugMessage("neurox::input::DataLoader::FilterLocalitySynapses...\n");
  CallAllNeurons(input::DataLoader::FilterRepeatedAndLinearizeContainers);
  DebugMessage("neurox::input::DataLoader::Finalize...\n");
  CallAllLocalities(input::DataLoader::Finalize);
  DebugMessage("neurox::Branch::BranchTree::InitLCOs...\n");
  CallAllNeurons(Branch::BranchTree::InitLCOs);
  input::Debugger::CompareMechanismsFunctions();
  input::Debugger::CompareAllBranches();

  if (neurox::input_params_->output_statistics_) {
    tools::Statistics::OutputMechanismsDistribution();
    tools::Statistics::OutputSimulationSize();
    // hpx_exit(0,NULL);
  }

  if (neurox::input_params_->load_balancing_)
    tools::LoadBalancing::PrintLoadBalancingTable();

  // call init action on each neuron (e.g. Finitialize, Cvodes init)
  DebugMessage("neurox::Branch::Initialize...\n");
  CallAllNeurons(Branch::Initialize);
#ifndef NDEBUG
  CallAllLocalities(input::Debugger::Finitialize);
  CallAllLocalities(input::Debugger::ThreadTableCheck);
  input::Debugger::CompareAllBranches();
#endif

  // iterator through all synchronizers (if many) and run
  const int synchronizer_id = (int)input_params_->synchronizer_;
  double tstop = input_params_->tstop_;

  printf("neurox::%s::%d neurons::starting...\n", synchronizer_->GetString(),
         neurox::neurons_count_);

  CallAllLocalities(Synchronizer::CallInitLocality, &synchronizer_id,
                    sizeof(synchronizer_id));
  CallAllNeurons(Synchronizer::CallInitNeuron);

  hpx_time_t time_now = hpx_time_now();
  if (input_params_->locality_comm_reduce_ || input_params_->scheduler_)
    CallAllLocalities(Synchronizer::RunLocality, &tstop, sizeof(tstop));
  else
    CallAllNeurons(Synchronizer::RunNeuron, &tstop, sizeof(tstop));
  double time_elapsed = hpx_time_elapsed_ms(time_now) / 1e3;

  printf(
      "neurox::%s: %d neurons, biological time: %.04f secs, solver time: "
      "%.02f secs\n",
      synchronizer_->GetString(), neurox::neurons_count_,
      input_params_->tstop_ / 1000., time_elapsed);

#ifdef NDEBUG
  // output benchmark info
  printf("csv,%d,%d,%d,%.1f,%.1f,%d,%d,%d,%.2f,%d,%.2f,%d,%.2f,%d,%.3f,%.3f\n",
         neurox::neurons_count_, hpx_get_num_ranks(), hpx_get_num_threads(),
         neurox::neurons_count_ / (double)hpx_get_num_ranks(),
         input_params_->tstop_, synchronizer_->GetId(),
         input_params_->graph_mechs_parallelism_ ? 1 : 0,
         input_params_->mech_instances_parallelism_ ? 1 : 0,
         input_params_->mech_instance_percent_per_block,
         input_params_->branch_parallelism_ ? 1 : 0,
         input_params_->subtree_complexity,
         input_params_->load_balancing_ ? 1 : 0,
         input_params_->subsection_complexity,
         input_params_->locality_comm_reduce_ ? 1 : 0,
         Mechanism::time_spent_in_mechs_ / 1e9,  // nano secs to secs
         time_elapsed);
  fflush(stdout);
#endif

  if (input_params_->output_comm_count_) {
    Statistics::CommCount::Counts counts;
    Statistics::CommCount::ReduceCounts(&counts);
    printf(
        "neurox::Statistics::CommCount:: p2p: %d; reduce:%d; %s; Avg spike "
        "rate: %.3f Hz/neuron\n",
        counts.point_to_point_count, counts.reduce_count,
        input_params_->locality_comm_reduce_ ? "comm-reduce" : "",
        (floble_t)counts.spike_count / (input_params_->tstop_ / 1000.) /
            (floble_t)neurox::neurons_count_);
  }

  CallAllNeurons(Synchronizer::CallClearNeuron);
  CallAllLocalities(Synchronizer::CallClearLocality);

  input::Debugger::RunCoreneuronAndCompareAllBranches();
  CallAllNeurons(Branch::Clear);
  hpx_bcast_rsync(neurox::Clear);

  double total_elapsed_time = hpx_time_elapsed_ms(total_time_now) / 1e3;
  DebugMessage(string("neurox::total time: " +
                      std::to_string(total_elapsed_time) + " secs\n")
                   .c_str());
  hpx_exit(0, NULL);
}

hpx_action_t Clear = 0;
int Clear_handler() {
  NEUROX_MEM_PIN(uint64_t);
  delete[] neurox::mechanisms_;
  delete[] neurox::neurons_;
  delete[] neurox::mechanisms_map_;
  delete synchronizer_;

  if (input_params_->locality_comm_reduce_ || input_params_->scheduler_) {
    (*neurox::locality::neurons_).clear();
    (*neurox::locality::netcons_branches_).clear();
    (*neurox::locality::netcons_somas_).clear();
    delete neurox::locality::neurons_;
    delete neurox::locality::netcons_branches_;
    delete neurox::locality::netcons_somas_;
    neurox::locality::netcons_branches_ = nullptr;
    neurox::locality::netcons_somas_ = nullptr;
    neurox::locality::neurons_ = nullptr;
  }

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

};  // namespace neurox
