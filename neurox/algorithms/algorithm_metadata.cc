#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::algorithms;

AlgorithmMetadata* AlgorithmMetadata::New(AlgorithmType type) {
  switch (type) {
    case AlgorithmType::kBackwardEulerDebug:
      return new DebugAlgorithm::CommunicationBarrier();
    case AlgorithmType::kBackwardEulerCoreneuron:
      return new CoreneuronAlgorithm::CommunicationBarrier();
    case AlgorithmType::kBackwardEulerAllReduce:
      return new AllreduceAlgorithm::AllReducesInfo();
    case AlgorithmType::kBackwardEulerSlidingTimeWindow:
      return new AllreduceAlgorithm::AllReducesInfo();
    case AlgorithmType::kBackwardEulerTimeDependencyLCO:
      return new TimeDependencyLCOAlgorithm::TimeDependencies();
    default:
      return nullptr;
  }
  return nullptr;
}
