#include "neurox/algorithms/sliding_time_window_algorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

hpx_t* SlidingTimeWindowAlgorithm::allReduces = nullptr;

SlidingTimeWindowAlgorithm::SlidingTimeWindowAlgorithm() {
  AllReduceAlgorithm::AllReducesInfo::reductionsPerCommStep =
      SlidingTimeWindowAlgorithm::kAllReducesCount;
}

SlidingTimeWindowAlgorithm::~SlidingTimeWindowAlgorithm() {}

const AlgorithmType SlidingTimeWindowAlgorithm::GetType() {
  return AlgorithmType::kBackwardEulerSlidingTimeWindow;
}

const char* SlidingTimeWindowAlgorithm::GetTypeString() {
  return "BackwardEulerSlidingTimeWindow";
}

void SlidingTimeWindowAlgorithm::Init() {
  AllReduceAlgorithm::SubscribeAllReduces(
      SlidingTimeWindowAlgorithm::allReduces,
      SlidingTimeWindowAlgorithm::kAllReducesCount);
}

void SlidingTimeWindowAlgorithm::Clear() {
  AllReduceAlgorithm::UnsubscribeAllReduces(
      SlidingTimeWindowAlgorithm::allReduces,
      SlidingTimeWindowAlgorithm::kAllReducesCount);
}

double SlidingTimeWindowAlgorithm::Launch() {
  int totalSteps = Algorithm::getTotalStepsCount();
  hpx_time_t now = hpx_time_now();
  if (input_params->all_reduce_at_locality)
    hpx_bcast_rsync(Branch::BackwardEulerOnLocality, &totalSteps, sizeof(int));
  else
    neurox::wrappers::CallAllNeurons(Branch::BackwardEuler, &totalSteps,
                                     sizeof(int));
  double elapsedTime = hpx_time_elapsed_ms(now) / 1e3;
  input::Debugger::RunCoreneuronAndCompareAllBranches();
  return elapsedTime;
}

void SlidingTimeWindowAlgorithm::StepBegin(Branch*) {}

void SlidingTimeWindowAlgorithm::StepEnd(Branch* b, hpx_t spikesLco) {
  AllReduceAlgorithm::WaitForSpikesDelivery(b, spikesLco);
  input::Debugger::SingleNeuronStepAndCompare(&nrn_threads[b->nt_->id], b,
                                              input_params->second_order_);
}

void SlidingTimeWindowAlgorithm::Run(Branch* b, const void* args) {
  AllReduceAlgorithm::Run2(b, args);
}

hpx_t SlidingTimeWindowAlgorithm::SendSpikes(Neuron* neuron, double tt,
                                             double) {
  return AllReduceAlgorithm::SendSpikes2(neuron, tt);
}
