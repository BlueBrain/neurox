#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::algorithms;

AlgorithmMetaData* AlgorithmMetaData::New(AlgorithmType type) {
  switch (type) {
    case AlgorithmType::BackwardEulerDebug:
      return new DebugAlgorithm::CommunicationBarrier();
    case AlgorithmType::BackwardEulerCoreneuron:
      return new CoreneuronAlgorithm::CommunicationBarrier();
    case AlgorithmType::BackwardEulerAllReduce:
      return new AllReduceAlgorithm::AllReducesInfo();
    case AlgorithmType::BackwardEulerSlidingTimeWindow:
      return new AllReduceAlgorithm::AllReducesInfo();
    case AlgorithmType::BackwardEulerTimeDependencyLCO:
      return new TimeDependencyLCOAlgorithm::TimeDependencies();
    default:
      return nullptr;
  }
  return nullptr;
}
