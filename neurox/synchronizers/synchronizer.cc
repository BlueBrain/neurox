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
    case Synchronizers::kTimeDependency:
      return new TimeDependencySynchronizer();
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
    case Synchronizers::kTimeDependency:
      return new TimeDependencySynchronizer::TimeDependencies();
    default:
      return nullptr;
  }
  return nullptr;
}
