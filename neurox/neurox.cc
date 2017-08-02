#include "neurox/neurox.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <map>

#include "coreneuron/nrniv/nrn_stats.h"
#include "coreneuron/nrniv/nrniv_decl.h"

using namespace neurox::algorithms;

namespace neurox {

hpx_t *neurons_ = nullptr;
int neurons_count_ = 0;
int mechanisms_count_ = -1;
int *mechanisms_map_ = nullptr;
neurox::Mechanism **mechanisms_ = nullptr;
neurox::tools::CmdLineParser *input_params_ = nullptr;
neurox::algorithms::Algorithm *algorithm_ = nullptr;

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

  if (neurox::input_params_->output_statistics_) {
    tools::Statistics::OutputMechanismsDistribution();
    tools::Statistics::OutputSimulationSize();
    // hpx_exit(0,NULL);
  }

  neurox::input::Debugger::CompareMechanismsFunctions();
  neurox::input::Debugger::CompareAllBranches();

  DebugMessage("neurox::Branch::Finitialize...\n");
  neurox::wrappers::CallAllNeurons(Branch::Finitialize);
#ifndef NDEBUG
  hpx_bcast_rsync(neurox::input::Debugger::Finitialize);
  neurox::input::Debugger::CompareAllBranches();
#endif

  DebugMessage("neurox::Branch::threadTableCheck...\n");
  neurox::wrappers::CallAllNeurons(Branch::ThreadTableCheck);
#ifndef NDEBUG
  hpx_bcast_rsync(neurox::input::Debugger::ThreadTableCheck);
  neurox::input::Debugger::CompareAllBranches();
#endif

  double total_time_elapsed = 0;
  if (input_params_->algorithm_ == AlgorithmType::kBenchmarkAll) {
    // TODO for this to work, we have to re-set algorothm in all cpus?
    for (int type = 0; type < 4; type++) {
      algorithm_ = Algorithm::New((AlgorithmType)type);
      algorithm_->Init();
      algorithm_->PrintStartInfo();
      double timeElapsed = algorithm_->Launch();
      total_time_elapsed += timeElapsed;
      algorithm_->Clear();
      delete algorithm_;

#ifdef NDEBUG
      // output benchmark info
      printf("csv,%d,%d,%d,%.1f,%.1f,%d,%d,%d,%d,%.2f\n", neurox::neurons_count,
             hpx_get_num_ranks(), hpx_get_num_threads(),
             neurox::neurons_count / (double)hpx_get_num_ranks(),
             input_params->tstop, algorithm->getType(),
             input_params->multiMex ? 1 : 0, input_params->branchingDepth,
             input_params->allReduceAtLocality ? 1 : 0, timeElapsed);
      fflush(stdout);
#endif
    }
  } else {
    algorithm_ = Algorithm::New(input_params_->algorithm_);
    algorithm_->Init();
    algorithm_->PrintStartInfo();
    total_time_elapsed = algorithm_->Launch();
    algorithm_->Clear();
    delete algorithm_;
  }

  printf(
      "neurox::end (%d neurons, biological time: %.3f secs, solver time: %.3f "
      "secs).\n",
      neurox::neurons_count_, input_params_->tstop_ / 1000.0,
      total_time_elapsed);

  neurox::wrappers::CallAllNeurons(Branch::Clear);
  hpx_bcast_rsync(neurox::Clear);
  hpx_exit(0, NULL);
}

hpx_action_t Clear = 0;
int Clear_handler() {
  NEUROX_MEM_PIN(uint64_t);
  delete[] neurox::mechanisms_;
  delete[] neurox::neurons_;
  delete[] neurox::mechanisms_map_;

  if (input_params_->all_reduce_at_locality_) {
    AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::localityNeurons
        ->clear();
    delete AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::
        localityNeurons;
    AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::localityNeurons =
        nullptr;
  }

#ifndef NDEBUG
  neurox::input::DataLoader::CleanCoreneuronData(true);
#endif
  return neurox::wrappers::MemoryUnpin(target);
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
