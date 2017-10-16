#include "neurox/neurox.h"

#include <algorithm>
#include <cstring>
#include <utility>

using namespace neurox;
using namespace neurox::synchronizers;
using namespace neurox::interpolators;

Neuron::Neuron(neuron_id_t neuron_id, floble_t ap_threshold)
    : gid_(neuron_id),
      threshold_(ap_threshold),
      synchronizer_neuron_info_(nullptr),
      synapses_localities_(nullptr) {
  this->synapses_transmission_flag_ = false;
  this->synapses_mutex_ = hpx_lco_sema_new(1);
  this->refractory_period_ = 0;

  /* Some synchronizers (eg Time Dependency) populate neuron
   * info during DataLoader, so must be instantiated now */
  SynchronizerIds id = (SynchronizerIds)input_params_->synchronizer_;
  if (id == SynchronizerIds::kBenchmarkAll)
    id = SynchronizerIds::kTimeDependency;  // First in Benchmark
  this->synchronizer_neuron_info_ = SynchronizerNeuronInfo::New(id);

  if (input_params_->locality_comm_reduce_)
    this->synapses_localities_ = new set<hpx_t>();
}

Neuron::~Neuron() {
  if (synapses_mutex_ != HPX_NULL) hpx_lco_delete_sync(synapses_mutex_);
  for (Synapse*& s : synapses_) delete s;
  delete synchronizer_neuron_info_;
}

Neuron::Synapse::Synapse(hpx_t branchAddr, floble_t minDelay,
                         hpx_t topBranchAddr, int destinationGid)
    : branch_addr_(branchAddr),
      min_delay_(minDelay),
      top_branch_addr_(topBranchAddr) {
  const double& teps = TimeDependencySynchronizer::TimeDependencies::kTEps;
  const double& notification_ratio =
      TimeDependencySynchronizer::TimeDependencies::kNotificationIntervalRatio;
  this->next_notification_time_ =
      input_params_->tstart_ + teps + this->min_delay_ * notification_ratio;
  this->previous_spike_lco_ = hpx_lco_future_new(0);
  hpx_lco_set_rsync(
      this->previous_spike_lco_, 0,
      NULL);  // starts as set and will be reset when synapses happen
#ifndef NDEBUG
  this->destination_gid_ = destinationGid;
#endif
}

Neuron::Synapse::~Synapse() {
  if (previous_spike_lco_ != HPX_NULL) {
    hpx_lco_delete_sync(previous_spike_lco_);
  }
}

size_t Neuron::GetSynapsesCount() {
  hpx_lco_sema_p(synapses_mutex_);
  size_t size = synapses_.size();
  if (synapses_localities_) {
    synapses_localities_->clear();
    synapses_localities_ = nullptr;
  }
  hpx_lco_sema_v_sync(synapses_mutex_);
  return size;
}

void Neuron::AddSynapse(Synapse* syn, hpx_t source_locality_addr) {
  hpx_lco_sema_p(synapses_mutex_);
  synapses_.push_back(syn);
  synapses_.shrink_to_fit();
  if (input_params_->locality_comm_reduce_)
    synapses_localities_->insert(source_locality_addr);
  hpx_lco_sema_v_sync(synapses_mutex_);
}

// netcvode.cpp::static bool pscheck(...)
bool Neuron::CheckAPthresholdAndTransmissionFlag(floble_t v) {
  // can only spike if AP threshold has been reach and spikes havent already
  // been transmitted
  if (v > threshold_) {
    if (synapses_transmission_flag_ == false) {
      synapses_transmission_flag_ = true;
      return true;
    }
  } else {
    synapses_transmission_flag_ = false;
  }
  return false;
}

hpx_t Neuron::SendSpikes(floble_t t)  // netcvode.cpp::PreSyn::send()
{
  const spike_time_t tt =
      (spike_time_t)t + 1e-10;  // Coreneuron logic, do not change!
#if !defined(NDEBUG)
  printf("== Neuron %d spiked at %.3f ms\n", this->gid_, tt);
#endif

  const bool locality_reduce = input_params_->locality_comm_reduce_;
  if ((!locality_reduce && synapses_.size() == 0) ||
      (locality_reduce && synapses_localities_->size() == 0))
    return HPX_NULL;
  return synchronizer_->SendSpikes(this, tt, t);
}
