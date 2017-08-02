#include "neurox/neurox.h"

#include <algorithm>
#include <cstring>
#include <utility>

using namespace neurox;
using namespace neurox::solver;
using namespace neurox::algorithms;

Neuron::Neuron(neuron_id_t neuron_id, floble_t ap_threshold)
    : gid_(neuron_id), threshold_(ap_threshold), algorithm_metadata_(nullptr) {
  this->synapses_transmission_flag_ = false;
  this->synapses_mutex_ = hpx_lco_sema_new(1);
  this->refractory_period_ = 0;
  this->algorithm_metadata_ = AlgorithmMetadata::New(input_params_->algorithm_);
  assert(this->algorithm_metadata_ != nullptr);
  assert(
      TimeDependencyLCOAlgorithm::TimeDependencies::notificationIntervalRatio >
          0 &&
      TimeDependencyLCOAlgorithm::TimeDependencies::notificationIntervalRatio <=
          1);
  assert(CoreneuronAlgorithm::CommunicationBarrier::kCommStepSize %
             AllReduceAlgorithm::AllReducesInfo::reductionsPerCommStep ==
         0);
}

Neuron::~Neuron() {
  if (synapses_mutex_ != HPX_NULL) hpx_lco_delete_sync(synapses_mutex_);
  for (Synapse*& s : synapses_) delete s;
  delete algorithm_metadata_;
}

Neuron::Synapse::Synapse(hpx_t branchAddr, floble_t minDelay,
                         hpx_t topBranchAddr, int destinationGid)
    : branch_addr_(branchAddr), min_delay_(minDelay), top_branch_addr_(topBranchAddr) {
  const double& teps = TimeDependencyLCOAlgorithm::TimeDependencies::teps;
  const double& notification_ratio =
      TimeDependencyLCOAlgorithm::TimeDependencies::notificationIntervalRatio;
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
  if (previous_spike_lco_ != HPX_NULL)
      hpx_lco_delete_sync(previous_spike_lco_);
}

size_t Neuron::GetSynapsesCount() {
  hpx_lco_sema_p(synapses_mutex_);
  size_t size = synapses_.size();
  hpx_lco_sema_v_sync(synapses_mutex_);
  return size;
}

void Neuron::AddSynapse(Synapse* syn) {
  hpx_lco_sema_p(synapses_mutex_);
  synapses_.push_back(syn);
  synapses_.shrink_to_fit();
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
  printf("== Neuron gid %d spiked at %.3f ms\n", this->gid_, tt);
#endif

  if (synapses_.size() == 0) return HPX_NULL;
  return algorithm_->SendSpikes(this, tt, t);
}
