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

void SlidingTimeWindowSynchronizer::NeuronReduceInit(Branch* b) {
  AllreduceSynchronizer::NeuronReduce(b, kAllReducesCount);
}

hpx_t SlidingTimeWindowSynchronizer::SendSpikes(Neuron* n, double tt, double) {
  return AllreduceSynchronizer::SendSpikes2(n, tt);
}

double SlidingTimeWindowSynchronizer::NeuronReduceInterval(Branch* b) {
  return AllreduceSynchronizer::NeuronReduceInterval2(b, kAllReducesCount);
}

double SlidingTimeWindowSynchronizer::LocalityReductionInterval() {
  return AllreduceSynchronizer::LocalityReduceInterval2(kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::LocalityReduceInit() {
  AllreduceSynchronizer::AllReduceLocalityInfo::LocalityReduce(
      kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::NeuronReduceEnd(Branch* b, hpx_t spikesLco) {
  AllreduceSynchronizer::WaitForSpikesDelivery(b, spikesLco);
}
