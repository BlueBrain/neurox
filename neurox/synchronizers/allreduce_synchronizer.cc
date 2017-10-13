#include "neurox/synchronizers/allreduce_synchronizer.h"

using namespace neurox;
using namespace neurox::synchronizers;
using namespace neurox::interpolators;

hpx_t* AllreduceSynchronizer::allreduces_ = nullptr;

AllreduceSynchronizer::AllreduceSynchronizer() {
  AllreduceSynchronizer::AllReducesInfo::reductions_per_comm_step_ =
      AllreduceSynchronizer::kAllReducesCount;
}

AllreduceSynchronizer::~AllreduceSynchronizer() {}

const Synchronizers AllreduceSynchronizer::GetId() {
  return Synchronizers::kAllReduce;
}

const char* AllreduceSynchronizer::GetString() {
  return "BackwardEulerAllReduce";
}

void AllreduceSynchronizer::Init() {
  SubscribeAllReduces(AllreduceSynchronizer::allreduces_,
                      AllreduceSynchronizer::kAllReducesCount);
}

void AllreduceSynchronizer::Clear() {
  UnsubscribeAllReduces(AllreduceSynchronizer::allreduces_,
                        AllreduceSynchronizer::kAllReducesCount);
}

void AllreduceSynchronizer::Launch() {
    /*
  if (input_params_->locality_comm_reduce_)
    hpx_bcast_rsync(BackwardEuler::RunOnLocality);
  else
    neurox::wrappers::CallAllNeurons(BackwardEuler::RunOnNeuron);
  input::Debugger::RunCoreneuronAndCompareAllBranches();
  */
}

void AllreduceSynchronizer::BeforeStep(Branch*) {}

void AllreduceSynchronizer::AfterStep(Branch* b, hpx_t spikesLco) {
  WaitForSpikesDelivery(b, spikesLco);
  input::Debugger::SingleNeuronStepAndCompare(&nrn_threads[b->nt_->id], b,
                                              input_params_->second_order_);
}

void AllreduceSynchronizer::Run(Branch* b, const void* args) { Run2(b, args); }

hpx_t AllreduceSynchronizer::SendSpikes(Neuron* n, double tt, double) {
  return Neuron::SendSpikesAsync(n, tt);
}

void AllreduceSynchronizer::SubscribeAllReduces(hpx_t*& allreduces,
                                                size_t allreduces_count) {
  assert(allreduces == nullptr);
  allreduces = new hpx_t[allreduces_count];

  hpx_bcast_rsync(
      AllreduceSynchronizer::AllReducesInfo::SetReductionsPerCommStep,
      &allreduces_count, sizeof(int));

  for (int i = 0; i < allreduces_count; i++)
    allreduces[i] = hpx_process_collective_allreduce_new(
        0, AllreduceSynchronizer::AllReducesInfo::Init,
        AllreduceSynchronizer::AllReducesInfo::Reduce);

  if (input_params_->locality_comm_reduce_)
    hpx_bcast_rsync(AllreduceSynchronizer::AllReducesInfo::AllReduceLocality::
                        SubscribeAllReduce,
                    allreduces, sizeof(hpx_t) * allreduces_count);
  else
    neurox::wrappers::CallAllNeurons(
        AllreduceSynchronizer::AllReducesInfo::SubscribeAllReduce, allreduces,
        sizeof(hpx_t) * allreduces_count);

  for (int i = 0; i < allreduces_count; i++)
    hpx_process_collective_allreduce_subscribe_finalize(allreduces[i]);
}

void AllreduceSynchronizer::UnsubscribeAllReduces(hpx_t*& allreduces,
                                                  size_t allreduces_count) {
  assert(allreduces != nullptr);
  if (input_params_->locality_comm_reduce_)
    hpx_bcast_rsync(AllreduceSynchronizer::AllReducesInfo::AllReduceLocality::
                        UnsubscribeAllReduce,
                    allreduces, sizeof(hpx_t) * allreduces_count);
  else
    neurox::wrappers::CallAllNeurons(
        AllreduceSynchronizer::AllReducesInfo::UnsubscribeAllReduce, allreduces,
        sizeof(hpx_t) * allreduces_count);

  for (int i = 0; i < allreduces_count; i++)
    hpx_process_collective_allreduce_delete(allreduces[i]);

  delete[] allreduces;
  allreduces = nullptr;
}

void AllreduceSynchronizer::WaitForSpikesDelivery(Branch* b, hpx_t spikes_lco) {
  // wait for spikes sent 4 steps ago (queue has always size 3)
  if (b->soma_) {
    AllReducesInfo* stw = (AllReducesInfo*)b->soma_->synchronizer_metadata_;
    assert(stw->spikes_lco_queue_.size() == BackwardEuler::GetMinSynapticDelaySteps() - 1);
    stw->spikes_lco_queue_.push(spikes_lco);
    hpx_t queued_spikes_lco = stw->spikes_lco_queue_.front();
    stw->spikes_lco_queue_.pop();
    if (queued_spikes_lco != HPX_NULL) {
      hpx_lco_wait(queued_spikes_lco);
      hpx_lco_delete_sync(queued_spikes_lco);
    }
  }
}

void AllreduceSynchronizer::Run2(Branch* b, const void* args) {
  int steps = *(int*)args;
  const int reductions_per_comm_step =
      AllreduceSynchronizer::AllReducesInfo::reductions_per_comm_step_;
  const int comm_step_size = BackwardEuler::GetMinSynapticDelaySteps();
  const int steps_per_reduction = comm_step_size / reductions_per_comm_step;
  const AllReducesInfo* stw =
      b->soma_ ? (AllReducesInfo*)b->soma_->synchronizer_metadata_ : nullptr;

  for (int s = 0; s < steps; s += comm_step_size)  // for every comm step
  {
#ifndef NDEBUG
    if (hpx_get_my_rank() == 0 && b->nt_->id == 0 && b->soma_) {
      printf("-- t=%.4f ms\n", input_params_->dt_ * s);
      fflush(stdout);
    }
#endif
    // for every reduction step
    for (int r = 0; r < reductions_per_comm_step; r++)  
    {
      if (b->soma_) {
        if (s >= comm_step_size)  // first comm-window does not wait
          hpx_lco_wait_reset(stw->allreduce_future_[r]);
        else
          // fixes crash for Synchronizer::ALL when running two hpx-reduce
          // -based
          // synchronizers in a row
          hpx_lco_reset_sync(stw->allreduce_future_[r]);

        hpx_process_collective_allreduce_join(stw->allreduce_lco_[r],
                                              stw->allreduce_id_[r], NULL, 0);
      }

      for (int n = 0; n < steps_per_reduction; n++)
          BackwardEuler::FullStep(b);

      // Input::Coreneuron::Debugger::stepAfterStepBackwardEuler(local,
      // &nrn_threads[this->nt->id], secondorder); //SMP ONLY
    }
  }
}

AllreduceSynchronizer::AllReducesInfo::AllReducesInfo() {
  for (int s = 0; s < BackwardEuler::GetMinSynapticDelaySteps() - 1; s++)
    this->spikes_lco_queue_.push(HPX_NULL);
}

AllreduceSynchronizer::AllReducesInfo::~AllReducesInfo() {
  for (int i = 0; i < spikes_lco_queue_.size(); i++) {
    hpx_t queued_spikes_lco = spikes_lco_queue_.front();
    if (queued_spikes_lco != HPX_NULL) hpx_lco_delete_sync(queued_spikes_lco);
    spikes_lco_queue_.pop();
  }
}

hpx_action_t AllreduceSynchronizer::AllReducesInfo::SubscribeAllReduce = 0;
int AllreduceSynchronizer::AllReducesInfo::SubscribeAllReduce_handler(
    const hpx_t* allreduces, const size_t size) {
  NEUROX_MEM_PIN(Branch);
  AllReducesInfo* stw = (AllReducesInfo*)local->soma_->synchronizer_metadata_;
  stw->allreduce_future_ = new hpx_t[AllReducesInfo::reductions_per_comm_step_];
  stw->allreduce_lco_ = new hpx_t[AllReducesInfo::reductions_per_comm_step_];
  stw->allreduce_id_ = new int[AllReducesInfo::reductions_per_comm_step_];
  for (int i = 0; i < size / sizeof(hpx_t); i++) {
    stw->allreduce_lco_[i] = allreduces[i];
    stw->allreduce_future_[i] =
        hpx_lco_future_new(0);  // no value to be reduced
    stw->allreduce_id_[i] = hpx_process_collective_allreduce_subscribe(
        allreduces[i], hpx_lco_set_action, stw->allreduce_future_[i]);
  }
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t AllreduceSynchronizer::AllReducesInfo::UnsubscribeAllReduce = 0;
int AllreduceSynchronizer::AllReducesInfo::UnsubscribeAllReduce_handler(
    const hpx_t* allreduces, const size_t size) {
  NEUROX_MEM_PIN(Branch);
  AllReducesInfo* stw = (AllReducesInfo*)local->soma_->synchronizer_metadata_;
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

int AllreduceSynchronizer::AllReducesInfo::reductions_per_comm_step_ = -1;
hpx_t* AllreduceSynchronizer::AllReducesInfo::AllReduceLocality::
    allreduce_future_ = nullptr;
hpx_t*
    AllreduceSynchronizer::AllReducesInfo::AllReduceLocality::allreduce_lco_ =
        nullptr;
int* AllreduceSynchronizer::AllReducesInfo::AllReduceLocality::allreduce_id_ =
    nullptr;

hpx_action_t AllreduceSynchronizer::AllReducesInfo::SetReductionsPerCommStep =
    0;
int AllreduceSynchronizer::AllReducesInfo::SetReductionsPerCommStep_handler(
    const int* val, const size_t) {
  NEUROX_MEM_PIN(uint64_t);
  reductions_per_comm_step_ = *val;
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t AllreduceSynchronizer::AllReducesInfo::AllReduceLocality::
    SubscribeAllReduce = 0;
int AllreduceSynchronizer::AllReducesInfo::AllReduceLocality::
    SubscribeAllReduce_handler(const hpx_t* allreduces, const size_t size) {
  NEUROX_MEM_PIN(uint64_t);
  assert(input_params_->locality_comm_reduce_);
  AllReduceLocality::allreduce_lco_ =
      new hpx_t[AllReducesInfo::reductions_per_comm_step_];
  AllReduceLocality::allreduce_future_ =
      new hpx_t[AllReducesInfo::reductions_per_comm_step_];
  AllReduceLocality::allreduce_id_ =
      new int[AllReducesInfo::reductions_per_comm_step_];
  for (int i = 0; i < size / sizeof(hpx_t); i++) {
    allreduce_lco_[i] = allreduces[i];
    allreduce_future_[i] = hpx_lco_future_new(0);  // no value to be reduced
    allreduce_id_[i] = hpx_process_collective_allreduce_subscribe(
        allreduces[i], hpx_lco_set_action, allreduce_future_[i]);
  }
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t AllreduceSynchronizer::AllReducesInfo::AllReduceLocality::
    UnsubscribeAllReduce = 0;
int AllreduceSynchronizer::AllReducesInfo::AllReduceLocality::
    UnsubscribeAllReduce_handler(const hpx_t* allreduces, const size_t size) {
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

hpx_action_t AllreduceSynchronizer::AllReducesInfo::Init = 0;
void AllreduceSynchronizer::AllReducesInfo::Init_handler(void*, const size_t) {}

hpx_action_t AllreduceSynchronizer::AllReducesInfo::Reduce = 0;
void AllreduceSynchronizer::AllReducesInfo::Reduce_handler(void*, const void*,
                                                           const size_t) {}

void AllreduceSynchronizer::AllReducesInfo::RegisterHpxActions() {
  wrappers::RegisterSingleVarAction<hpx_t>(SubscribeAllReduce,
                                           SubscribeAllReduce_handler);
  wrappers::RegisterSingleVarAction<hpx_t>(UnsubscribeAllReduce,
                                           UnsubscribeAllReduce_handler);
  wrappers::RegisterSingleVarAction<hpx_t>(
      AllReduceLocality::SubscribeAllReduce,
      AllReduceLocality::SubscribeAllReduce_handler);
  wrappers::RegisterSingleVarAction<hpx_t>(
      AllReduceLocality::UnsubscribeAllReduce,
      AllReduceLocality::UnsubscribeAllReduce_handler);

  wrappers::RegisterSingleVarAction<int>(
      AllReducesInfo::SetReductionsPerCommStep,
      AllReducesInfo::SetReductionsPerCommStep_handler);
  wrappers::RegisterAllReduceInitAction<void>(AllReducesInfo::Init,
                                              AllReducesInfo::Init_handler);
  wrappers::RegisterAllReduceReduceAction<void>(AllReducesInfo::Reduce,
                                                AllReducesInfo::Reduce_handler);
}
