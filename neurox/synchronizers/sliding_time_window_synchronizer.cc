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
  AllreduceSynchronizer::SubscribeAllReduces(
      SlidingTimeWindowSynchronizer::kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::Clear() {
  AllreduceSynchronizer::UnsubscribeAllReduces(
      SlidingTimeWindowSynchronizer::kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::BeforeStep(Branch*) {}

void SlidingTimeWindowSynchronizer::AfterStep(Branch* b, hpx_t spikesLco) {
  AllreduceSynchronizer::WaitForSpikesDelivery(b, spikesLco);
  input::Debugger::SingleNeuronStepAndCompare(&nrn_threads[b->nt_->id], b,
                                              input_params_->second_order_);
}

void SlidingTimeWindowSynchronizer::Run(Branch* b, const void* args) {
  AllreduceSynchronizer::Run2(b, args);
}

hpx_t SlidingTimeWindowSynchronizer::SendSpikes(Neuron* n, double tt, double) {
  return Neuron::SendSpikesAsync(n, tt);
}

double SlidingTimeWindowSynchronizer::GetLocalityReductionInterval()
{
    return neurox::min_synaptic_delay_/kAllReducesCount;
}

void SlidingTimeWindowSynchronizer::LocalityReduce() {
    AllreduceSynchronizer::AllReduceLocalityInfo::LocalityReduce(
                SlidingTimeWindowSynchronizer::kAllReducesCount
                );
}
