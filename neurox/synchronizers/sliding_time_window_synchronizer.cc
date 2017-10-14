#include "neurox/synchronizers/sliding_time_window_synchronizer.h"

using namespace neurox;
using namespace neurox::synchronizers;

SlidingTimeWindowSynchronizer::SlidingTimeWindowSynchronizer() {}

SlidingTimeWindowSynchronizer::~SlidingTimeWindowSynchronizer() {}

const SynchronizerIds SlidingTimeWindowSynchronizer::GetId() {
  return SynchronizerIds::kSlidingTimeWindow;
}

const char* SlidingTimeWindowSynchronizer::GetString() {
  return "BackwardEulerSlidingTimeWindow";
}

void SlidingTimeWindowSynchronizer::InitLocality() {
  AllreduceSynchronizer::SubscribeAllReduces(kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::ClearLocality() {
  AllreduceSynchronizer::UnsubscribeAllReduces(kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::BeforeSteps(Branch* branch) {
  AllreduceSynchronizer::NeuronReduce(branch, kAllReducesCount);
}

double SlidingTimeWindowSynchronizer::GetMaxStepTime(Branch* b) {
  AllreduceSynchronizer::GetMaxStepTime2(b, kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::AfterSteps(Branch* b, hpx_t spikesLco) {
  AllreduceSynchronizer::WaitForSpikesDelivery(b, spikesLco);
  input::Debugger::SingleNeuronStepAndCompare(&nrn_threads[b->nt_->id], b,
                                              input_params_->second_order_);
}

hpx_t SlidingTimeWindowSynchronizer::SendSpikes(Neuron* n, double tt, double) {
  return Neuron::SendSpikesAsync(n, tt);
}

double SlidingTimeWindowSynchronizer::GetLocalityReductionInterval() {
  return AllreduceSynchronizer::GetLocalityReductionInterval2(kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::LocalityReduce() {
  AllreduceSynchronizer::AllReduceLocalityInfo::LocalityReduce(
      kAllReducesCount);
}
