#include "neurox/neurox.h"

#include <algorithm>
#include <cstring>
#include <utility>

using namespace neurox;
using namespace neurox::synchronizers;
using namespace neurox::interpolators;

Neuron::Neuron(neuron_id_t neuron_id, floble_t ap_threshold)
    : containers_buffer_(nullptr),
      containers_buffer_size_(0),
      gid_(neuron_id),
      threshold_(ap_threshold),
      synapses_linear_(nullptr),
      synchronizer_neuron_info_(nullptr),
      scheduler_step_trigger_(HPX_NULL) {
  this->synapses_transmission_flag_ = false;
  this->synapses_mutex_ = hpx_lco_sema_new(1);
  this->refractory_period_ = 0;
  SynchronizerIds id = (SynchronizerIds)input_params_->synchronizer_;
  this->synchronizer_neuron_info_ = SynchronizerNeuronInfo::New(id);
}

Neuron::~Neuron() {
  for (Synapse*& s : synapses_) delete s;
  delete synchronizer_neuron_info_;
  if (synapses_mutex_ != HPX_NULL) hpx_lco_delete_sync(synapses_mutex_);
  if (scheduler_step_trigger_ != HPX_NULL)
    hpx_lco_delete_sync(scheduler_step_trigger_);

  if (containers_buffer_) delete containers_buffer_;
}

Neuron::Synapse::Synapse(hpx_t branch_addr, floble_t min_delay,
                         hpx_t soma_or_locality_addr, int destination_gid)
    : branch_addr_(branch_addr),
      min_delay_(min_delay),
      soma_or_locality_addr_(soma_or_locality_addr) {
  const double& teps = TimeDependencySynchronizer::TimeDependencies::kTEps;
  const double& notification_ratio =
      TimeDependencySynchronizer::TimeDependencies::kNotificationIntervalRatio;
  this->next_notification_time_ =
      input_params_->tstart_ + teps + this->min_delay_ * notification_ratio;
  this->previous_synapse_lco_ = hpx_lco_future_new(0);
  // starts LCOs as set and will be reset when synapses happen
  hpx_lco_set_rsync(this->previous_synapse_lco_, 0, NULL);
  //#ifndef NDEBUG
  this->destination_gid_ = destination_gid;
  //#endif
}

Neuron::Synapse::~Synapse() {
  if (previous_synapse_lco_ != HPX_NULL)
    hpx_lco_delete_sync(previous_synapse_lco_);
}

size_t Neuron::GetSynapsesCount() {
  hpx_lco_sema_p(synapses_mutex_);
  size_t size = synapses_linear_ ? synapses_linear_->Count() : synapses_.size();
  hpx_lco_sema_v_sync(synapses_mutex_);
  return size;
}

Neuron::Synapse* Neuron::GetSynapseAtOffset(size_t i) {
  return synapses_linear_ ? synapses_linear_->At(i) : synapses_.at(i);
}

void Neuron::AddSynapse(Synapse* syn) {
  /* for locality-based reduction, repeated synapses
   * will be filtered by DataLoader::Finalize */
  hpx_lco_sema_p(synapses_mutex_);
  synapses_.push_back(syn);
  synapses_.shrink_to_fit();
  hpx_lco_sema_v_sync(synapses_mutex_);
}

void Neuron::LinearizeContainers() {
  assert(input_params_->synchronizer_ == SynchronizerIds::kTimeDependency);
  TimeDependencySynchronizer::TimeDependencies* td =
      (TimeDependencySynchronizer::TimeDependencies*)this
          ->synchronizer_neuron_info_;

  assert(td->dependencies_min_delay_.size() ==
         td->dependencies_max_time_allowed_.size());

  containers_buffer_size_ = 0;

  // in-place constructors
  containers_buffer_size_ +=
      Vectorizer::SizeOf(linear::Vector<Synapse>::Size(synapses_.size()));
  containers_buffer_size_ +=
      Vectorizer::SizeOf(linear::Map<neuron_id_t, floble_t>::Size(
          td->dependencies_min_delay_.size()));
  containers_buffer_size_ +=
      Vectorizer::SizeOf(linear::Map<neuron_id_t, floble_t>::Size(
          td->dependencies_max_time_allowed_.size()));
  containers_buffer_ = new unsigned char[containers_buffer_size_];

  size_t containers_buffer_it = 0;
  synapses_linear_ =
      (linear::Vector<Synapse>*)&containers_buffer_[containers_buffer_it];
  new (synapses_linear_)
      linear::Vector<Synapse>(synapses_, (unsigned char*)synapses_linear_);
  containers_buffer_it +=
      Vectorizer::SizeOf(linear::Vector<Synapse>::Size(synapses_.size()));

  td->dependencies_min_delay_linear_ =
      (linear::Map<neuron_id_t,
                   floble_t>*)&containers_buffer_[containers_buffer_it];
  new (td->dependencies_min_delay_linear_) linear::Map<neuron_id_t, floble_t>(
      td->dependencies_min_delay_,
      (unsigned char*)td->dependencies_min_delay_linear_);
  containers_buffer_it +=
      Vectorizer::SizeOf(linear::Map<neuron_id_t, floble_t>::Size(
          td->dependencies_min_delay_.size()));

  td->dependencies_max_time_allowed_linear_ =
      (linear::Map<neuron_id_t,
                   floble_t>*)&containers_buffer_[containers_buffer_it];
  new (td->dependencies_max_time_allowed_linear_)
      linear::Map<neuron_id_t, floble_t>(
          td->dependencies_max_time_allowed_,
          (unsigned char*)td->dependencies_max_time_allowed_linear_);
  containers_buffer_it +=
      Vectorizer::SizeOf(linear::Map<neuron_id_t, floble_t>::Size(
          td->dependencies_max_time_allowed_.size()));
  assert(containers_buffer_it == containers_buffer_size_);

#ifndef NDEBUG
  assert(synapses_.size() == synapses_linear_->Count());
  size_t syn_count = GetSynapsesCount();
  for (int i = 0; i < syn_count; i++) {
    Neuron::Synapse* s = GetSynapseAtOffset(i);
    Neuron::Synapse* s1 = synapses_.at(i);
    Neuron::Synapse* s2 = synapses_linear_->At(i);
    assert(s1->branch_addr_ == s2->branch_addr_);
    assert(s1->min_delay_ == s2->min_delay_);
    assert(s1->previous_synapse_lco_ == s2->previous_synapse_lco_);
    assert(s1->soma_or_locality_addr_ == s2->soma_or_locality_addr_);
  }
#endif
  synapses_.clear();
  td->dependencies_max_time_allowed_.clear();
  td->dependencies_min_delay_.clear();
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
#if !defined(NDEBUG) or defined(PRINT_NEURON_SPIKED)
  fprintf(stderr, "-- Neuron %d spiked at %.3f ms\n", this->gid_, tt);
#endif

  if (GetSynapsesCount() > 0) synchronizer_->SendSpikes(this, tt, t);
}
