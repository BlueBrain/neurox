/*
# =============================================================================
# Copyright (c) 2015 - 2021-2022 Blue Brain Project/EPFL
#
# See top-level LICENSE file for details.
# =============================================================================.
*/
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

void SlidingTimeWindowSynchronizer::NeuronSyncInit(Branch* b) {
  AllreduceSynchronizer::NeuronReduce(b, kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::SendSpikes(Neuron* n, double tt, double) {
  return AllreduceSynchronizer::SendSpikes2(n, tt);
}

double SlidingTimeWindowSynchronizer::GetNeuronMaxStep(Branch*) {
  return AllreduceSynchronizer::NeuronReduceInterval2(kAllReducesCount);
}

double SlidingTimeWindowSynchronizer::LocalitySyncInterval() {
  return AllreduceSynchronizer::LocalityReduceInterval2(kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::LocalitySyncInit() {
  AllreduceSynchronizer::AllReduceLocalityInfo::LocalityReduce(
      kAllReducesCount);
}

void SlidingTimeWindowSynchronizer::NeuronSyncEnd(Branch* b) {
  AllreduceSynchronizer::WaitForSpikesDelivery(b);
}
