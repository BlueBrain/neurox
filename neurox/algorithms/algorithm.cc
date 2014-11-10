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
      return new AllReduceAlgorithm();
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
         GetTypeString(), neurox::neurons_count, input_params->tstop / 1000,
         input_params->dt);
  fflush(stdout);
}

int Algorithm::getTotalStepsCount() {
  return (input_params->tstop - input_params->tstart) / input_params->dt;
}
