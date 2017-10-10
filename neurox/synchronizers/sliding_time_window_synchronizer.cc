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
  Synchronizer::FixedStepMethodsInit();
  AllreduceSynchronizer::SubscribeAllReduces(
      SlidingTimeWindowSynchronizer::allreduces_,
      SlidingTimeWindowSynchronizer::kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::Clear() {
  AllreduceSynchronizer::UnsubscribeAllReduces(
      SlidingTimeWindowSynchronizer::allreduces_,
      SlidingTimeWindowSynchronizer::kAllReducesCount);
}

double SlidingTimeWindowSynchronizer::Launch() {
  int total_steps = Synchronizer::GetTotalStepsCount();
  hpx_time_t now = hpx_time_now();
  if (input_params_->allreduce_at_locality_)
    hpx_bcast_rsync(Branch::BackwardEulerOnLocality, &total_steps, sizeof(int));
  else
    neurox::wrappers::CallAllNeurons(Branch::BackwardEuler, &total_steps,
                                     sizeof(int));
  double elapsed_time = hpx_time_elapsed_ms(now) / 1e3;
  input::Debugger::RunCoreneuronAndCompareAllBranches();
  return elapsed_time;
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
