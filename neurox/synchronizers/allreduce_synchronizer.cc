/*
# =============================================================================
# Copyright (c) 2015 - 2021-2022 Blue Brain Project/EPFL
#
# See top-level LICENSE file for details.
# =============================================================================.
*/
#include "neurox/synchronizers/allreduce_synchronizer.h"

using namespace neurox;
using namespace neurox::synchronizers;
using namespace neurox::interpolators;

hpx_t* AllreduceSynchronizer::allreduces_ = nullptr;
int AllreduceSynchronizer::AllReduceLocalityInfo::next_allreduce_id_ = -1;

AllreduceSynchronizer::AllreduceSynchronizer() {}

AllreduceSynchronizer::~AllreduceSynchronizer() {}

const SynchronizerIds AllreduceSynchronizer::GetId() {
  return SynchronizerIds::kAllReduce;
}

const char* AllreduceSynchronizer::GetString() {
  return "AllreduceSynchronizer";
}

void AllreduceSynchronizer::InitLocality() {
  SubscribeAllReduces(kAllReducesCount);
}

void AllreduceSynchronizer::ClearLocality() {
  UnsubscribeAllReduces(kAllReducesCount);
}

void AllreduceSynchronizer::NeuronSyncInit(Branch* b) {
  NeuronReduce(b, kAllReducesCount);
}

void AllreduceSynchronizer::SendSpikes(Neuron* n, double tt, double) {
  return SendSpikes2(n, tt);
}

double AllreduceSynchronizer::GetNeuronMaxStep(Branch*) {
  return AllreduceSynchronizer::NeuronReduceInterval2(kAllReducesCount);
}

double AllreduceSynchronizer::LocalitySyncInterval() {
  return AllreduceSynchronizer::LocalityReduceInterval2(kAllReducesCount);
}

void AllreduceSynchronizer::LocalitySyncInit() {
  AllReduceLocalityInfo::LocalityReduce(kAllReducesCount);
}

void AllreduceSynchronizer::NeuronSyncEnd(Branch* b) {
  WaitForSpikesDelivery(b);
}

void AllreduceSynchronizer::SendSpikes2(Neuron* neuron, spike_time_t tt) {
  hpx_action_t spike_action = input_params_->locality_comm_reduce_
                                  ? Branch::AddSpikeEventLocality
                                  : Branch::AddSpikeEvent;
  size_t syn_count = neuron->GetSynapsesCount();

  hpx_t new_synapses_lco = hpx_lco_and_new(syn_count);
  for (int i = 0; i < syn_count; i++) {
    Neuron::Synapse* s = neuron->GetSynapseAtOffset(i);
    hpx_call(s->branch_addr_, spike_action, new_synapses_lco, &neuron->gid_,
             sizeof(neuron_id_t), &tt, sizeof(spike_time_t));
  }

  AllReduceNeuronInfo* stw =
      (AllReduceNeuronInfo*)neuron->synchronizer_neuron_info_;
  stw->spikes_lco_queue_.push(
      std::make_pair(tt + neurox::min_synaptic_delay_, new_synapses_lco));
}

void AllreduceSynchronizer::NeuronReduce(const Branch* branch,
                                         const int allreduces_count) {
  // if locality-reduction is on, neurons no not participate n reduction
  if (input_params_->locality_comm_reduce_) return;
  assert(branch->soma_);

  if (input_params_->output_comm_count_) {
    hpx_lco_sema_p(Statistics::CommCount::mutex);
    Statistics::CommCount::counts.reduce_count++;
    hpx_lco_sema_v_sync(Statistics::CommCount::mutex);
  }

  AllReduceNeuronInfo* stw =
      (AllReduceNeuronInfo*)branch->soma_->synchronizer_neuron_info_;

  // if reduction id < 0, it's still on the first comm-window
  // (within first comm-window, synchronizer does not wait)
  int& r = stw->next_allreduce_id_;
  int allreduce_id = r < 0 ? r + allreduces_count : r;
  if (r < 0)
    /* fixes crash for Synchronizer::ALL when running two
     * hpx-reduce -based synchronizers in a row*/
    hpx_lco_reset_sync(stw->allreduce_future_[allreduce_id]);
  else
    /* neuron-level reduction */
    hpx_lco_wait_reset(stw->allreduce_future_[allreduce_id]);

  /* mark neuron's step in this all reduce */
  hpx_process_collective_allreduce_join(stw->allreduce_lco_[allreduce_id],
                                        stw->allreduce_id_[allreduce_id], NULL,
                                        0);

  if (++r == allreduces_count) r = 0;
}

double AllreduceSynchronizer::NeuronReduceInterval2(
    const int allreduces_count) {
  return neurox::min_synaptic_delay_ / allreduces_count;
}

double AllreduceSynchronizer::LocalityReduceInterval2(
    const double allreduces_count) {
  return neurox::min_synaptic_delay_ / allreduces_count;
}

void AllreduceSynchronizer::SubscribeAllReduces(size_t allreduces_count) {
  // rank 0 ask all ranks or neurons to unsubscribe
  if (hpx_get_my_rank() > 0) return;

  assert(allreduces_ == nullptr);
  allreduces_ = new hpx_t[allreduces_count];

  for (int i = 0; i < allreduces_count; i++)
    allreduces_[i] = hpx_process_collective_allreduce_new(
        0, AllReduceNeuronInfo::Init, AllReduceNeuronInfo::Reduce);

  if (input_params_->locality_comm_reduce_)
    hpx_bcast_rsync(AllReduceLocalityInfo::Subscribe, allreduces_,
                    sizeof(hpx_t) * allreduces_count);
  else
    wrappers::CallAllNeurons(AllReduceNeuronInfo::Subscribe, allreduces_,
                             sizeof(hpx_t) * allreduces_count);

  for (int i = 0; i < allreduces_count; i++)
    hpx_process_collective_allreduce_subscribe_finalize(allreduces_[i]);
}

void AllreduceSynchronizer::UnsubscribeAllReduces(size_t allreduces_count) {
  // rank 0 asks all ranks or neurons to unsubscribe
  if (hpx_get_my_rank() > 0) return;

  if (input_params_->locality_comm_reduce_)
    hpx_bcast_rsync(AllReduceLocalityInfo::Unsubscribe, allreduces_,
                    sizeof(hpx_t) * allreduces_count);
  else
    wrappers::CallAllNeurons(AllReduceNeuronInfo::Unsubscribe, allreduces_,
                             sizeof(hpx_t) * allreduces_count);

  for (int i = 0; i < allreduces_count; i++)
    hpx_process_collective_allreduce_delete(allreduces_[i]);

  delete[] allreduces_;
  allreduces_ = nullptr;
}

void AllreduceSynchronizer::WaitForSpikesDelivery(Branch* b) {
  // wait for spikes sent 4 steps ago (queue has always size 3)
  if (b->soma_) {
    const double t = b->nt_->_t;
    AllReduceNeuronInfo* stw =
        (AllReduceNeuronInfo*)b->soma_->synchronizer_neuron_info_;
    while (!stw->spikes_lco_queue_.empty() &&
           t + 0.01 > stw->spikes_lco_queue_.top().first) {
      TimedSpike queued_spike = stw->spikes_lco_queue_.top();
      stw->spikes_lco_queue_.pop();
      hpx_lco_wait(queued_spike.second);
      hpx_lco_delete_sync(queued_spike.second);
    }
  }
}

void AllreduceSynchronizer::AllReduceLocalityInfo::LocalityReduce(
    int allreduces_count) {
  // if no reduction at locality level, reduction is done by neurons
  if (!input_params_->locality_comm_reduce_) return;

  if (input_params_->output_comm_count_) {
    hpx_lco_sema_p(Statistics::CommCount::mutex);
    Statistics::CommCount::counts.reduce_count++;
    hpx_lco_sema_v_sync(Statistics::CommCount::mutex);
  }

  // if reduction id < 0, it's still on the first comm-window
  // within first comm-window, synchronizer does not wait
  int& r = next_allreduce_id_;
  int allreduce_id = r < 0 ? r + allreduces_count : r;
  if (r < 0)
    /* fixes crash for Synchronizer::ALL when running two
     *  hpx-reduce -based synchronizers in a row*/
    hpx_lco_reset_sync(allreduce_future_[allreduce_id]);
  else
    /* locality-level reduction */
    hpx_lco_wait_reset(allreduce_future_[allreduce_id]);

  /* mark locality step on this allreduce */
  hpx_process_collective_allreduce_join(allreduce_lco_[allreduce_id],
                                        allreduce_id_[allreduce_id], NULL, 0);

  if (++r == allreduces_count) r = 0;
}

AllreduceSynchronizer::AllReduceNeuronInfo::AllReduceNeuronInfo(
    const size_t allreduces_count) {
  // negative value means ignore until it reaches 0
  next_allreduce_id_ = -allreduces_count;
}

AllreduceSynchronizer::AllReduceNeuronInfo::~AllReduceNeuronInfo() {
  for (int i = 0; i < spikes_lco_queue_.size(); i++)
    hpx_lco_delete_sync(spikes_lco_queue_.top().second);
}

hpx_action_t AllreduceSynchronizer::AllReduceNeuronInfo::Subscribe = 0;
int AllreduceSynchronizer::AllReduceNeuronInfo::Subscribe_handler(
    const hpx_t* allreduces, const size_t size) {
  NEUROX_MEM_PIN(Branch);
  const int allreduces_count = size / sizeof(hpx_t);
  AllReduceNeuronInfo* stw =
      (AllReduceNeuronInfo*)local->soma_->synchronizer_neuron_info_;
  stw->allreduce_future_ = new hpx_t[allreduces_count];
  stw->allreduce_lco_ = new hpx_t[allreduces_count];
  stw->allreduce_id_ = new int[allreduces_count];
  for (int i = 0; i < allreduces_count; i++) {
    stw->allreduce_lco_[i] = allreduces[i];
    stw->allreduce_future_[i] =
        hpx_lco_future_new(0);  // no value to be reduced
    stw->allreduce_id_[i] = hpx_process_collective_allreduce_subscribe(
        allreduces[i], hpx_lco_set_action, stw->allreduce_future_[i]);
  }
  NEUROX_MEM_UNPIN;
}

hpx_action_t AllreduceSynchronizer::AllReduceNeuronInfo::Unsubscribe = 0;
int AllreduceSynchronizer::AllReduceNeuronInfo::Unsubscribe_handler(
    const hpx_t* allreduces, const size_t size) {
  NEUROX_MEM_PIN(Branch);
  AllReduceNeuronInfo* stw =
      (AllReduceNeuronInfo*)local->soma_->synchronizer_neuron_info_;
  for (int i = 0; i < size / sizeof(hpx_t); i++) {
    hpx_process_collective_allreduce_unsubscribe(allreduces[i],
                                                 stw->allreduce_id_[i]);
    if (stw->allreduce_future_[i] != HPX_NULL)
      hpx_lco_delete_sync(stw->allreduce_future_[i]);
  }
  delete[] stw->allreduce_lco_;
  stw->allreduce_lco_ = nullptr;
  delete[] stw->allreduce_future_;
  stw->allreduce_future_ = nullptr;
  delete[] stw->allreduce_id_;
  stw->allreduce_id_ = nullptr;
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_t* AllreduceSynchronizer::AllReduceLocalityInfo::allreduce_future_ =
    nullptr;
hpx_t* AllreduceSynchronizer::AllReduceLocalityInfo::allreduce_lco_ = nullptr;
int* AllreduceSynchronizer::AllReduceLocalityInfo::allreduce_id_ = nullptr;

hpx_action_t AllreduceSynchronizer::AllReduceLocalityInfo::Subscribe = 0;
int AllreduceSynchronizer::AllReduceLocalityInfo::Subscribe_handler(
    const hpx_t* allreduces, const size_t size) {
  NEUROX_MEM_PIN(uint64_t);
  assert(input_params_->locality_comm_reduce_);
  const int allreduces_count = size / sizeof(hpx_t);
  AllReduceLocalityInfo::allreduce_lco_ = new hpx_t[allreduces_count];
  AllReduceLocalityInfo::allreduce_future_ = new hpx_t[allreduces_count];
  AllReduceLocalityInfo::allreduce_id_ = new int[allreduces_count];
  for (int i = 0; i < allreduces_count; i++) {
    allreduce_lco_[i] = allreduces[i];
    allreduce_future_[i] = hpx_lco_future_new(0);  // no value to be reduced
    allreduce_id_[i] = hpx_process_collective_allreduce_subscribe(
        allreduces[i], hpx_lco_set_action, allreduce_future_[i]);
  }

  // negative value means ignore until it reaches 0
  AllReduceLocalityInfo::next_allreduce_id_ = -allreduces_count;
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t AllreduceSynchronizer::AllReduceLocalityInfo::Unsubscribe = 0;
int AllreduceSynchronizer::AllReduceLocalityInfo::Unsubscribe_handler(
    const hpx_t* allreduces, const size_t size) {
  NEUROX_MEM_PIN(uint64_t);
  assert(input_params_->locality_comm_reduce_);
  for (int i = 0; i < size / sizeof(hpx_t); i++) {
    hpx_process_collective_allreduce_unsubscribe(allreduces[i],
                                                 allreduce_id_[i]);
    hpx_lco_delete_sync(allreduce_future_[i]);
  }
  delete[] allreduce_lco_;
  allreduce_lco_ = nullptr;
  delete[] allreduce_future_;
  allreduce_future_ = nullptr;
  delete[] allreduce_id_;
  allreduce_id_ = nullptr;
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t AllreduceSynchronizer::AllReduceNeuronInfo::Init = 0;
void AllreduceSynchronizer::AllReduceNeuronInfo::Init_handler(void*,
                                                              const size_t) {}

hpx_action_t AllreduceSynchronizer::AllReduceNeuronInfo::Reduce = 0;
void AllreduceSynchronizer::AllReduceNeuronInfo::Reduce_handler(void*,
                                                                const void*,
                                                                const size_t) {}

void AllreduceSynchronizer::RegisterHpxActions() {
  // node level functions
  wrappers::RegisterSingleVarAction<hpx_t>(
      AllReduceNeuronInfo::Subscribe, AllReduceNeuronInfo::Subscribe_handler);
  wrappers::RegisterSingleVarAction<hpx_t>(
      AllReduceNeuronInfo::Unsubscribe,
      AllReduceNeuronInfo::Unsubscribe_handler);
  wrappers::RegisterAllReduceInitAction<void>(
      AllReduceNeuronInfo::Init, AllReduceNeuronInfo::Init_handler);
  wrappers::RegisterAllReduceReduceAction<void>(
      AllReduceNeuronInfo::Reduce, AllReduceNeuronInfo::Reduce_handler);

  // locality leve functions
  wrappers::RegisterSingleVarAction<hpx_t>(
      AllReduceLocalityInfo::Subscribe,
      AllReduceLocalityInfo::Subscribe_handler);
  wrappers::RegisterSingleVarAction<hpx_t>(
      AllReduceLocalityInfo::Unsubscribe,
      AllReduceLocalityInfo::Unsubscribe_handler);
}
