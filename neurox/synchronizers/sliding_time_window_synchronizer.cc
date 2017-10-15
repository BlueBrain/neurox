#include "neurox/synchronizers/sliding_time_window_synchronizer.h"

using namespace neurox;
using namespace neurox::synchronizers;

SlidingTimeWindowSynchronizer::SlidingTimeWindowSynchronizer() {}

SlidingTimeWindowSynchronizer::~SlidingTimeWindowSynchronizer() {}

const SynchronizerIds SlidingTimeWindowSynchronizer::GetId() {
  return SynchronizerIds::kSlidingTimeWindow;
}

const char* SlidingTimeWindowSynchronizer::GetString() {
  return "SlidingTimeWindowSynchronizer";
}

void SlidingTimeWindowSynchronizer::InitLocality() {
  AllreduceSynchronizer::SubscribeAllReduces(kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::ClearLocality() {
  AllreduceSynchronizer::UnsubscribeAllReduces(kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::BeforeSteps(Branch* b) {
 AllreduceSynchronizer:: NeuronReduce(b, kAllReducesCount);
}

double SlidingTimeWindowSynchronizer::GetMaxStep(Branch* b) {
  AllreduceSynchronizer::GetMaxStep2(b, kAllReducesCount);
}

double SlidingTimeWindowSynchronizer::GetLocalityReductionInterval() {
  return AllreduceSynchronizer::GetLocalityReductionInterval2(kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::LocalityReduce() {
  AllreduceSynchronizer::AllReduceLocalityInfo::LocalityReduce(kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::AfterSteps(Branch* b, hpx_t spikesLco) {
  AllreduceSynchronizer::WaitForSpikesDelivery(b, spikesLco);
}

hpx_t SlidingTimeWindowSynchronizer::SendSpikes(Neuron* n, double tt, double) {
  return Neuron::SendSpikesAsync(n, tt);
}
