#include "neurox/algorithms/sliding_time_window_algorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

hpx_t* SlidingTimeWindowAlgorithm::allreduces_ = nullptr;

SlidingTimeWindowAlgorithm::SlidingTimeWindowAlgorithm() {
  AllreduceAlgorithm::AllReducesInfo::reductions_per_comm_step_ =
      SlidingTimeWindowAlgorithm::kAllReducesCount;
}

SlidingTimeWindowAlgorithm::~SlidingTimeWindowAlgorithm() {}

const AlgorithmId SlidingTimeWindowAlgorithm::GetId() {
  return AlgorithmId::kBackwardEulerSlidingTimeWindow;
}

const char* SlidingTimeWindowAlgorithm::GetString() {
  return "BackwardEulerSlidingTimeWindow";
}

void SlidingTimeWindowAlgorithm::Init() {
    Algorithm::FixedStepMethodsInit();
  AllreduceAlgorithm::SubscribeAllReduces(
      SlidingTimeWindowAlgorithm::allreduces_,
      SlidingTimeWindowAlgorithm::kAllReducesCount);
}

void SlidingTimeWindowAlgorithm::Clear() {
  AllreduceAlgorithm::UnsubscribeAllReduces(
      SlidingTimeWindowAlgorithm::allreduces_,
      SlidingTimeWindowAlgorithm::kAllReducesCount);
}

double SlidingTimeWindowAlgorithm::Launch() {
  int total_steps = Algorithm::GetTotalStepsCount();
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

void SlidingTimeWindowAlgorithm::StepBegin(Branch*) {}

void SlidingTimeWindowAlgorithm::StepEnd(Branch* b, hpx_t spikesLco) {
  AllreduceAlgorithm::WaitForSpikesDelivery(b, spikesLco);
  input::Debugger::SingleNeuronStepAndCompare(&nrn_threads[b->nt_->id], b,
                                              input_params_->second_order_);
}

void SlidingTimeWindowAlgorithm::Run(Branch* b, const void* args) {
  AllreduceAlgorithm::Run2(b, args);
}

hpx_t SlidingTimeWindowAlgorithm::SendSpikes(Neuron* n, double tt,
                                             double) {
  return Neuron::SendSpikesAsync(n, tt);
}
