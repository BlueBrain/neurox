#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::synchronizers;

Synchronizer* Synchronizer::New(Synchronizers type) {
  switch (type) {
    case Synchronizers::kDebug:
      return new DebugSynchronizer();
    case Synchronizers::kCoreneuron:
      return new CoreneuronSynchronizer();
    case Synchronizers::kAllReduce:
      return new AllreduceSynchronizer();
    case Synchronizers::kSlidingTimeWindow:
      return new SlidingTimeWindowSynchronizer();
    case Synchronizers::kTimeDependencyLCO:
      return new TimeDependencyLCOSynchronizer();
    default:
      return nullptr;
  }
  return nullptr;
};

SynchronizerMetadata* SynchronizerMetadata::New(Synchronizers type) {
  switch (type) {
    case Synchronizers::kDebug:
      return new DebugSynchronizer::CommunicationBarrier();
    case Synchronizers::kCoreneuron:
      return new CoreneuronSynchronizer::CommunicationBarrier();
    case Synchronizers::kAllReduce:
      return new AllreduceSynchronizer::AllReducesInfo();
    case Synchronizers::kSlidingTimeWindow:
      return new AllreduceSynchronizer::AllReducesInfo();
    case Synchronizers::kTimeDependencyLCO:
      return new TimeDependencyLCOSynchronizer::TimeDependencies();
    default:
      return nullptr;
  }
  return nullptr;
}

void Synchronizer::PrintStartInfo() {
  printf("neurox::Synchronizer::%s (%d neurons, t=%.03f secs, dt=%.03f milisecs\n",
         GetString(), neurox::neurons_count_, input_params_->tstop_ / 1000,
         input_params_->dt_);
  fflush(stdout);
}

int Synchronizer::GetTotalStepsCount() {
  return (input_params_->tstop_ - input_params_->tstart_) / input_params_->dt_;
}

void Synchronizer::FixedStepMethodsInit() {
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
