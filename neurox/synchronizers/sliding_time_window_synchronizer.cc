#include "neurox/synchronizers/sliding_time_window_synchronizer.h"

using namespace neurox;
using namespace neurox::synchronizers;

hpx_t* SlidingTimeWindowSynchronizer::allreduces_ = nullptr;

SlidingTimeWindowSynchronizer::SlidingTimeWindowSynchronizer() {
  AllreduceSynchronizer::AllReducesInfo::reductions_per_comm_step_ =
      SlidingTimeWindowSynchronizer::kAllReducesCount;
}

SlidingTimeWindowSynchronizer::~SlidingTimeWindowSynchronizer() {}

const Synchronizers SlidingTimeWindowSynchronizer::GetId() {
  return Synchronizers::kSlidingTimeWindow;
}

const char* SlidingTimeWindowSynchronizer::GetString() {
  return "BackwardEulerSlidingTimeWindow";
}

void SlidingTimeWindowSynchronizer::Init() {
  AllreduceSynchronizer::SubscribeAllReduces(
      SlidingTimeWindowSynchronizer::allreduces_,
      SlidingTimeWindowSynchronizer::kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::Clear() {
  AllreduceSynchronizer::UnsubscribeAllReduces(
      SlidingTimeWindowSynchronizer::allreduces_,
      SlidingTimeWindowSynchronizer::kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::Launch() {
    /*
  int total_steps = interpolators::BackwardEuler::GetTotalStepsCount();
  if (input_params_->locality_comm_reduce_)
    hpx_bcast_rsync(interpolators::BackwardEuler::RunOnLocality, &total_steps, sizeof(int));
  else
    neurox::wrappers::CallAllNeurons(interpolators::BackwardEuler::RunOnNeuron, &total_steps,
                                     sizeof(int));
  input::Debugger::RunCoreneuronAndCompareAllBranches();
  */
}

void SlidingTimeWindowSynchronizer::StepBegin(Branch*) {}

void SlidingTimeWindowSynchronizer::StepEnd(Branch* b, hpx_t spikesLco) {
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
