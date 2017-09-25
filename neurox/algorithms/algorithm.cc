#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::algorithms;

Algorithm* Algorithm::New(AlgorithmId type) {
  switch (type) {
    case AlgorithmId::kBackwardEulerDebug:
      return new DebugAlgorithm();
    case AlgorithmId::kBackwardEulerCoreneuron:
      return new CoreneuronAlgorithm();
    case AlgorithmId::kBackwardEulerAllReduce:
      return new AllreduceAlgorithm();
    case AlgorithmId::kBackwardEulerSlidingTimeWindow:
      return new SlidingTimeWindowAlgorithm();
    case AlgorithmId::kBackwardEulerTimeDependencyLCO:
      return new TimeDependencyLCOAlgorithm();
    case AlgorithmId::kCvodes:
      return new CvodesAlgorithm();
    default:
      return nullptr;
  }
  return nullptr;
};

AlgorithmMetadata* AlgorithmMetadata::New(AlgorithmId type) {
  switch (type) {
    case AlgorithmId::kBackwardEulerDebug:
      return new DebugAlgorithm::CommunicationBarrier();
    case AlgorithmId::kBackwardEulerCoreneuron:
      return new CoreneuronAlgorithm::CommunicationBarrier();
    case AlgorithmId::kBackwardEulerAllReduce:
      return new AllreduceAlgorithm::AllReducesInfo();
    case AlgorithmId::kBackwardEulerSlidingTimeWindow:
      return new AllreduceAlgorithm::AllReducesInfo();
    case AlgorithmId::kBackwardEulerTimeDependencyLCO:
      return new TimeDependencyLCOAlgorithm::TimeDependencies();
    case AlgorithmId::kCvodes:
      return new CvodesAlgorithm::BranchCvodes();
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
