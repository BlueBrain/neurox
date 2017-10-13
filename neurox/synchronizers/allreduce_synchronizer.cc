#include "neurox/synchronizers/allreduce_synchronizer.h"

using namespace neurox;
using namespace neurox::synchronizers;
using namespace neurox::interpolators;

hpx_t* AllreduceSynchronizer::allreduces_ = nullptr;
int AllreduceSynchronizer::AllReduceLocalityInfo::next_allreduce_id=-1;

AllreduceSynchronizer::AllreduceSynchronizer() {}

AllreduceSynchronizer::~AllreduceSynchronizer() {}

const Synchronizers AllreduceSynchronizer::GetId() {
  return Synchronizers::kAllReduce;
}

const char* AllreduceSynchronizer::GetString() {
  return "BackwardEulerAllReduce";
}

void AllreduceSynchronizer::Init() {
  SubscribeAllReduces(AllreduceSynchronizer::kAllReducesCount);
}

void AllreduceSynchronizer::Clear() {
  UnsubscribeAllReduces(AllreduceSynchronizer::kAllReducesCount);
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

double AllreduceSynchronizer::GetLocalityReductionInterval()
{
    return neurox::min_synaptic_delay_/kAllReducesCount;
}

void AllreduceSynchronizer::LocalityReduce() {
    AllreduceSynchronizer::AllReduceLocalityInfo::LocalityReduce(
                AllreduceSynchronizer::kAllReducesCount
                );
}


void AllreduceSynchronizer::SubscribeAllReduces(size_t allreduces_count) {
  assert(allreduces_ == nullptr);
  allreduces_ = new hpx_t[allreduces_count];

  for (int i = 0; i < allreduces_count; i++)
    allreduces_[i] = hpx_process_collective_allreduce_new(
        0, AllReducesInfo::Init, AllReducesInfo::Reduce);

  if (input_params_->locality_comm_reduce_)
    hpx_bcast_rsync(AllReduceLocalityInfo::SubscribeAllReduce,
                    allreduces_, sizeof(hpx_t) * allreduces_count);
  else
    neurox::wrappers::CallAllNeurons(
        AllReducesInfo::SubscribeAllReduce, allreduces_,
        sizeof(hpx_t) * allreduces_count);

  for (int i = 0; i < allreduces_count; i++)
    hpx_process_collective_allreduce_subscribe_finalize(allreduces_[i]);
}

void AllreduceSynchronizer::UnsubscribeAllReduces(size_t allreduces_count) {
  assert(allreduces_ != nullptr);
  if (input_params_->locality_comm_reduce_)
    hpx_bcast_rsync(AllReduceLocalityInfo::UnsubscribeAllReduce,
                    allreduces_, sizeof(hpx_t) * allreduces_count);
  else
    neurox::wrappers::CallAllNeurons(
        AllReducesInfo::UnsubscribeAllReduce, allreduces_,
        sizeof(hpx_t) * allreduces_count);

  for (int i = 0; i < allreduces_count; i++)
    hpx_process_collective_allreduce_delete(allreduces_[i]);

  delete[] allreduces_;
  allreduces_ = nullptr;
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
  const int reductions_per_comm_step = -1;
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

void AllreduceSynchronizer::AllReduceLocalityInfo::LocalityReduce(int allreduces_count) {

    // if reduction id < 0, it's still on the first comm-window
    //within first comm-window, synchronizer does not wait
    int r = next_allreduce_id;
    if (r>=0)
      hpx_lco_wait_reset(allreduce_future_[r]);
    else
      // fixes crash for Synchronizer::ALL when running two hpx-reduce -based
      // synchronizers in a row
      hpx_lco_reset_sync(allreduce_future_[r]);

    hpx_process_collective_allreduce_join(
                allreduce_lco_[r], allreduce_id_[r], NULL, 0);

    if (++next_allreduce_id == allreduces_count)
        next_allreduce_id=0;
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
  const int allreduces_count = size / sizeof(hpx_t);
  AllReducesInfo* stw = (AllReducesInfo*)local->soma_->synchronizer_metadata_;
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

hpx_t* AllreduceSynchronizer::AllReduceLocalityInfo::
    allreduce_future_ = nullptr;
hpx_t*
    AllreduceSynchronizer::AllReduceLocalityInfo::allreduce_lco_ =
        nullptr;
int* AllreduceSynchronizer::AllReduceLocalityInfo::allreduce_id_ =
    nullptr;

hpx_action_t AllreduceSynchronizer::AllReduceLocalityInfo::
    SubscribeAllReduce = 0;
int AllreduceSynchronizer::AllReduceLocalityInfo::
    SubscribeAllReduce_handler(const hpx_t* allreduces, const size_t size) {
  NEUROX_MEM_PIN(uint64_t);
  assert(input_params_->locality_comm_reduce_);
  const int allreduces_count = size / sizeof(hpx_t);
  AllReduceLocalityInfo::allreduce_lco_ =
      new hpx_t[allreduces_count];
  AllReduceLocalityInfo::allreduce_future_ =
      new hpx_t[allreduces_count];
  AllReduceLocalityInfo::allreduce_id_ =
      new int[allreduces_count];
  for (int i = 0; i < allreduces_count; i++) {
    allreduce_lco_[i] = allreduces[i];
    allreduce_future_[i] = hpx_lco_future_new(0);  // no value to be reduced
    allreduce_id_[i] = hpx_process_collective_allreduce_subscribe(
        allreduces[i], hpx_lco_set_action, allreduce_future_[i]);
  }

  //negative value means ignore until it reaches 0
  AllReduceLocalityInfo::next_allreduce_id = -allreduces_count;
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t AllreduceSynchronizer::AllReduceLocalityInfo::
    UnsubscribeAllReduce = 0;
int AllreduceSynchronizer::AllReduceLocalityInfo::
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
      AllReduceLocalityInfo::SubscribeAllReduce,
      AllReduceLocalityInfo::SubscribeAllReduce_handler);
  wrappers::RegisterSingleVarAction<hpx_t>(
      AllReduceLocalityInfo::UnsubscribeAllReduce,
      AllReduceLocalityInfo::UnsubscribeAllReduce_handler);

  wrappers::RegisterAllReduceInitAction<void>(AllReducesInfo::Init,
                                              AllReducesInfo::Init_handler);
  wrappers::RegisterAllReduceReduceAction<void>(AllReducesInfo::Reduce,
                                                AllReducesInfo::Reduce_handler);
}
