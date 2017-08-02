#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::algorithms;

Algorithm* Algorithm::New(AlgorithmType type) {
  switch (type) {
    case AlgorithmType::kBackwardEulerDebug:
      return new DebugAlgorithm();
    case AlgorithmType::kBackwardEulerCoreneuron:
      return new CoreneuronAlgorithm();
    case AlgorithmType::kBackwardEulerAllReduce:
      return new AllreduceAlgorithm();
    case AlgorithmType::kBackwardEulerSlidingTimeWindow:
      return new SlidingTimeWindowAlgorithm();
    case AlgorithmType::kBackwardEulerTimeDependencyLCO:
      return new TimeDependencyLCOAlgorithm();
    default:
      return nullptr;
  }
  return nullptr;
};

void Algorithm::PrintStartInfo() {
  printf("neurox::Algorithm::%s (%d neurons, t=%.03f secs, dt=%.03f milisecs\n",
         GetTypeString(), neurox::neurons_count_, input_params_->tstop_ / 1000,
         input_params_->dt_);
  fflush(stdout);
}

int Algorithm::GetTotalStepsCount() {
  return (input_params_->tstop_ - input_params_->tstart_) / input_params_->dt_;
}
