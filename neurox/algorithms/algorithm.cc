#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::algorithms;

Algorithm* Algorithm::New(SyncAlgorithms type) {
  switch (type) {
    case SyncAlgorithms::kDebug:
      return new DebugAlgorithm();
    case SyncAlgorithms::kCoreneuronMPICollective:
      return new CoreneuronAlgorithm();
    case SyncAlgorithms::kAllReduce:
      return new AllreduceAlgorithm();
    case SyncAlgorithms::kSlidingTimeWindow:
      return new SlidingTimeWindowAlgorithm();
    case SyncAlgorithms::kTimeDependencyLCO:
      return new TimeDependencyLCOAlgorithm();
    default:
      return nullptr;
  }
  return nullptr;
};

AlgorithmMetadata* AlgorithmMetadata::New(SyncAlgorithms type) {
  switch (type) {
    case SyncAlgorithms::kDebug:
      return new DebugAlgorithm::CommunicationBarrier();
    case SyncAlgorithms::kCoreneuronMPICollective:
      return new CoreneuronAlgorithm::CommunicationBarrier();
    case SyncAlgorithms::kAllReduce:
      return new AllreduceAlgorithm::AllReducesInfo();
    case SyncAlgorithms::kSlidingTimeWindow:
      return new AllreduceAlgorithm::AllReducesInfo();
    case SyncAlgorithms::kTimeDependencyLCO:
      return new TimeDependencyLCOAlgorithm::TimeDependencies();
    default:
      return nullptr;
  }
  return nullptr;
}

void Algorithm::PrintStartInfo() {
  printf("neurox::Algorithm::%s (%d neurons, t=%.03f secs, dt=%.03f milisecs\n",
         GetString(), neurox::neurons_count_, input_params_->tstop_ / 1000,
         input_params_->dt_);
  fflush(stdout);
}

int Algorithm::GetTotalStepsCount() {
  return (input_params_->tstop_ - input_params_->tstart_) / input_params_->dt_;
}

void Algorithm::FixedStepMethodsInit() {
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
}
