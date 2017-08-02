#include "neurox/neurox.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <map>

#include "coreneuron/nrniv/nrn_stats.h"
#include "coreneuron/nrniv/nrniv_decl.h"

using namespace neurox::algorithms;

namespace neurox {

/*
template<typename T, typename... ArgTypes>
int CountVars(T t, ArgTypes... args)
{
  return 1+CountVars(args...);
}
template<> int CountVars() {return 0;}

template<typename... Args>
hpx_status_t CallAllNeurons(hpx_action_t f, Args... args)
{
  va_list ap;
  hpx_t lco= hpx_lco_and_new(neurox::neurons_count);
  int e = HPX_SUCCESS;
  int n = neurox::CountVars(args...);
  for (size_t i = 0; i < neurox::neurons_count; i++)
    e += _hpx_call(neurox::neurons[i], f, lco, n, args...);
  hpx_lco_wait_reset(lco);
  hpx_lco_delete_sync(lco);
  return e;
}
*/

hpx_t *neurons = nullptr;
int neurons_count = 0;
int mechanisms_count = -1;
int *mechanisms_map = nullptr;
neurox::Mechanism **mechanisms = nullptr;
neurox::tools::CmdLineParser *input_params = nullptr;
neurox::algorithms::Algorithm *algorithm = nullptr;

Mechanism *GetMechanismFromType(int type) {
  assert(mechanisms_map[type] != -1);
  return mechanisms[mechanisms_map[type]];
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

  if (neurox::input_params->output_statistics_) {
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

  double totalTimeElapsed = 0;
  if (input_params->algorithm == AlgorithmType::kBenchmarkAll) {
    // TODO for this to work, we have to re-set algorothm in all cpus?
    for (int type = 0; type < 4; type++) {
      algorithm = Algorithm::New((AlgorithmType)type);
      algorithm->Init();
      algorithm->PrintStartInfo();
      double timeElapsed = algorithm->Launch();
      totalTimeElapsed += timeElapsed;
      algorithm->Clear();
      delete algorithm;

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
    algorithm = Algorithm::New(input_params->algorithm);
    algorithm->Init();
    algorithm->PrintStartInfo();
    totalTimeElapsed = algorithm->Launch();
    algorithm->Clear();
    delete algorithm;
  }

  printf(
      "neurox::end (%d neurons, biological time: %.3f secs, solver time: %.3f "
      "secs).\n",
      neurox::neurons_count, input_params->tstop_ / 1000.0, totalTimeElapsed);

  neurox::wrappers::CallAllNeurons(Branch::Clear);
  hpx_bcast_rsync(neurox::Clear);
  hpx_exit(0, NULL);
}

hpx_action_t Clear = 0;
int Clear_handler() {
  NEUROX_MEM_PIN(uint64_t);
  delete[] neurox::mechanisms;
  delete[] neurox::neurons;
  delete[] neurox::mechanisms_map;

  if (input_params->all_reduce_at_locality) {
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
