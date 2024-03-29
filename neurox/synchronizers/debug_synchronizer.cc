/*
# =============================================================================
# Copyright (c) 2015 - 2021-2022 Blue Brain Project/EPFL
#
# See top-level LICENSE file for details.
# =============================================================================.
*/
#include "neurox/synchronizers/debug_synchronizer.h"

using namespace neurox;
using namespace neurox::synchronizers;
using namespace neurox::interpolators;

DebugSynchronizer::DebugSynchronizer() {}

DebugSynchronizer::~DebugSynchronizer() {}

DebugSynchronizer::CommunicationBarrier::CommunicationBarrier() {
  this->all_spikes_lco_ = HPX_NULL;
}

DebugSynchronizer::CommunicationBarrier::~CommunicationBarrier() {
  if (all_spikes_lco_ != HPX_NULL) hpx_lco_delete_sync(all_spikes_lco_);
}

const SynchronizerIds DebugSynchronizer::GetId() {
  return SynchronizerIds::kDebug;
}

const char* DebugSynchronizer::GetString() { return "DebugSynchronizer"; }

void DebugSynchronizer::InitLocality() {}

void DebugSynchronizer::ClearLocality() {}

double DebugSynchronizer::LocalitySyncInterval() {
  return neurox::min_synaptic_delay_;
}

double DebugSynchronizer::GetNeuronMaxStep(Branch* b) { return b->nt_->_dt; }

void DebugSynchronizer::NeuronSyncInit(Branch*) {}

void DebugSynchronizer::NeuronSyncEnd(Branch* b) {
  if (b->soma_)  // end of comm-step (steps is the number of steps per commSize)
  {
    CommunicationBarrier* comm_barrier =
        (CommunicationBarrier*)b->soma_->synchronizer_neuron_info_;
    if (comm_barrier->all_spikes_lco_ != HPX_NULL)  // was set/used once
      hpx_lco_wait(comm_barrier->all_spikes_lco_);  // wait if needed
  }
}

void DebugSynchronizer::SendSpikes(Neuron* neuron, double tt, double) {
  CommunicationBarrier* comm_barrier =
      (CommunicationBarrier*)neuron->synchronizer_neuron_info_;

  hpx_t& all_spikes_lco = comm_barrier->all_spikes_lco_;

  if (all_spikes_lco == HPX_NULL)  // first use
    all_spikes_lco = hpx_lco_and_new(neuron->synapses_.size());
  else
    hpx_lco_reset_sync(all_spikes_lco);  // reset to use after

  AllreduceSynchronizer::SendSpikes2(neuron, tt);
}
