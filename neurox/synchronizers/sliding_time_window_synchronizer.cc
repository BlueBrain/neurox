#include "neurox/synchronizers/sliding_time_window_synchronizer.h"

using namespace neurox;
using namespace neurox::synchronizers;

SlidingTimeWindowSynchronizer::SlidingTimeWindowSynchronizer() {}

SlidingTimeWindowSynchronizer::~SlidingTimeWindowSynchronizer() {}

const Synchronizers SlidingTimeWindowSynchronizer::GetId() {
  return Synchronizers::kSlidingTimeWindow;
}

const char* SlidingTimeWindowSynchronizer::GetString() {
  return "BackwardEulerSlidingTimeWindow";
}

void SlidingTimeWindowSynchronizer::Init() {
  AllreduceSynchronizer::SubscribeAllReduces(kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::Clear() {
  AllreduceSynchronizer::UnsubscribeAllReduces(kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::BeforeStep(Branch* branch) {
  AllreduceSynchronizer::BeforeStep2(branch, kAllReducesCount);
}

double SlidingTimeWindowSynchronizer::GetMaxStepTime(Branch* b) {
  AllreduceSynchronizer::GetMaxStepTime2(b, kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::AfterStep(Branch* b, hpx_t spikesLco) {
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
