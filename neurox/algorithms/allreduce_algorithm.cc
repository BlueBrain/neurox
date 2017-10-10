#include "neurox/algorithms/allreduce_algorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

hpx_t* AllreduceAlgorithm::allreduces_ = nullptr;

AllreduceAlgorithm::AllreduceAlgorithm() {
  AllreduceAlgorithm::AllReducesInfo::reductions_per_comm_step_ =
      AllreduceAlgorithm::kAllReducesCount;
}

AllreduceAlgorithm::~AllreduceAlgorithm() {}

const Algorithms AllreduceAlgorithm::GetId() { return Algorithms::kAllReduce; }

const char* AllreduceAlgorithm::GetString() { return "BackwardEulerAllReduce"; }

void AllreduceAlgorithm::Init() {
  Algorithm::FixedStepMethodsInit();
  SubscribeAllReduces(AllreduceAlgorithm::allreduces_,
                      AllreduceAlgorithm::kAllReducesCount);
}

void AllreduceAlgorithm::Clear() {
  UnsubscribeAllReduces(AllreduceAlgorithm::allreduces_,
                        AllreduceAlgorithm::kAllReducesCount);
}

double AllreduceAlgorithm::Launch() {
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

void AllreduceAlgorithm::StepBegin(Branch*) {}

void AllreduceAlgorithm::StepEnd(Branch* b, hpx_t spikesLco) {
  WaitForSpikesDelivery(b, spikesLco);
  input::Debugger::SingleNeuronStepAndCompare(&nrn_threads[b->nt_->id], b,
                                              input_params_->second_order_);
}

void AllreduceAlgorithm::Run(Branch* b, const void* args) { Run2(b, args); }

hpx_t AllreduceAlgorithm::SendSpikes(Neuron* n, double tt, double) {
  return Neuron::SendSpikesAsync(n, tt);
}

void AllreduceAlgorithm::SubscribeAllReduces(hpx_t*& allreduces,
                                             size_t allreduces_count) {
  assert(allreduces == nullptr);
  allreduces = new hpx_t[allreduces_count];

  hpx_bcast_rsync(AllreduceAlgorithm::AllReducesInfo::SetReductionsPerCommStep,
                  &allreduces_count, sizeof(int));

  for (int i = 0; i < allreduces_count; i++)
    allreduces[i] = hpx_process_collective_allreduce_new(
        0, AllreduceAlgorithm::AllReducesInfo::Init,
        AllreduceAlgorithm::AllReducesInfo::Reduce);

  if (input_params_->allreduce_at_locality_)
    hpx_bcast_rsync(AllreduceAlgorithm::AllReducesInfo::AllReduceLocality::
                        SubscribeAllReduce,
                    allreduces, sizeof(hpx_t) * allreduces_count);
  else
    neurox::wrappers::CallAllNeurons(
        AllreduceAlgorithm::AllReducesInfo::SubscribeAllReduce, allreduces,
        sizeof(hpx_t) * allreduces_count);

  for (int i = 0; i < allreduces_count; i++)
    hpx_process_collective_allreduce_subscribe_finalize(allreduces[i]);
}

void AllreduceAlgorithm::UnsubscribeAllReduces(hpx_t*& allreduces,
                                               size_t allreduces_count) {
  assert(allreduces != nullptr);
  if (input_params_->allreduce_at_locality_)
    hpx_bcast_rsync(AllreduceAlgorithm::AllReducesInfo::AllReduceLocality::
                        UnsubscribeAllReduce,
                    allreduces, sizeof(hpx_t) * allreduces_count);
  else
    neurox::wrappers::CallAllNeurons(
        AllreduceAlgorithm::AllReducesInfo::UnsubscribeAllReduce, allreduces,
        sizeof(hpx_t) * allreduces_count);

  for (int i = 0; i < allreduces_count; i++)
    hpx_process_collective_allreduce_delete(allreduces[i]);

  delete[] allreduces;
  allreduces = nullptr;
}

void AllreduceAlgorithm::WaitForSpikesDelivery(Branch* b, hpx_t spikes_lco) {
  // wait for spikes sent 4 steps ago (queue has always size 3)
  if (b->soma_) {
    AllReducesInfo* stw = (AllReducesInfo*)b->soma_->algorithm_metadata_;
    assert(stw->spikes_lco_queue_.size() == neurox::min_delay_steps_ - 1);
    stw->spikes_lco_queue_.push(spikes_lco);
    hpx_t queued_spikes_lco = stw->spikes_lco_queue_.front();
    stw->spikes_lco_queue_.pop();
    if (queued_spikes_lco != HPX_NULL) {
      hpx_lco_wait(queued_spikes_lco);
      hpx_lco_delete_sync(queued_spikes_lco);
    }
  }
}

void AllreduceAlgorithm::Run2(Branch* b, const void* args) {
  int steps = *(int*)args;
  const int reductions_per_comm_step =
      AllreduceAlgorithm::AllReducesInfo::reductions_per_comm_step_;
  const int comm_step_size = neurox::min_delay_steps_;
  const int steps_per_reduction = comm_step_size / reductions_per_comm_step;
  const AllReducesInfo* stw =
      b->soma_ ? (AllReducesInfo*)b->soma_->algorithm_metadata_ : nullptr;

  for (int s = 0; s < steps;
       s += comm_step_size)  // for every communication step
  {
#ifndef NDEBUG
    if (hpx_get_my_rank() == 0 && b->nt_->id == 0 && b->soma_) {
      printf("-- t=%.4f ms\n", input_params_->dt_ * s);
      fflush(stdout);
    }
#endif
    for (int r = 0; r < reductions_per_comm_step;
         r++)  // for every reduction step
    {
      if (b->soma_) {
        if (s >= comm_step_size)  // first comm-window does not wait
          hpx_lco_wait_reset(stw->allreduce_future_[r]);
        else
          // fixes crash for Algorithm::ALL when running two hpx-reduce -based
          // algorithms in a row
          hpx_lco_reset_sync(stw->allreduce_future_[r]);

        hpx_process_collective_allreduce_join(stw->allreduce_lco_[r],
                                              stw->allreduce_id_[r], NULL, 0);
      }

      for (int n = 0; n < steps_per_reduction; n++) b->BackwardEulerStep();
      // Input::Coreneuron::Debugger::stepAfterStepBackwardEuler(local,
      // &nrn_threads[this->nt->id], secondorder); //SMP ONLY
    }
  }
}

AllreduceAlgorithm::AllReducesInfo::AllReducesInfo() {
  for (int s = 0; s < neurox::min_delay_steps_ - 1; s++)
    this->spikes_lco_queue_.push(HPX_NULL);
}

AllreduceAlgorithm::AllReducesInfo::~AllReducesInfo() {
  for (int i = 0; i < spikes_lco_queue_.size(); i++) {
    hpx_t queued_spikes_lco = spikes_lco_queue_.front();
    if (queued_spikes_lco != HPX_NULL) hpx_lco_delete_sync(queued_spikes_lco);
    spikes_lco_queue_.pop();
  }
}

hpx_action_t AllreduceAlgorithm::AllReducesInfo::SubscribeAllReduce = 0;
int AllreduceAlgorithm::AllReducesInfo::SubscribeAllReduce_handler(
    const hpx_t* allreduces, const size_t size) {
  NEUROX_MEM_PIN(Branch);
  AllReducesInfo* stw = (AllReducesInfo*)local->soma_->algorithm_metadata_;
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

hpx_action_t AllreduceAlgorithm::AllReducesInfo::UnsubscribeAllReduce = 0;
int AllreduceAlgorithm::AllReducesInfo::UnsubscribeAllReduce_handler(
    const hpx_t* allreduces, const size_t size) {
  NEUROX_MEM_PIN(Branch);
  AllReducesInfo* stw = (AllReducesInfo*)local->soma_->algorithm_metadata_;
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

int AllreduceAlgorithm::AllReducesInfo::reductions_per_comm_step_ = -1;
std::vector<hpx_t>*
    AllreduceAlgorithm::AllReducesInfo::AllReduceLocality::locality_neurons_ =
        nullptr;
hpx_t*
    AllreduceAlgorithm::AllReducesInfo::AllReduceLocality::allreduce_future_ =
        nullptr;
hpx_t* AllreduceAlgorithm::AllReducesInfo::AllReduceLocality::allreduce_lco_ =
    nullptr;
int* AllreduceAlgorithm::AllReducesInfo::AllReduceLocality::allreduce_id_ =
    nullptr;

hpx_action_t AllreduceAlgorithm::AllReducesInfo::SetReductionsPerCommStep = 0;
int AllreduceAlgorithm::AllReducesInfo::SetReductionsPerCommStep_handler(
    const int* val, const size_t) {
  NEUROX_MEM_PIN(uint64_t);
  reductions_per_comm_step_ = *val;
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t
    AllreduceAlgorithm::AllReducesInfo::AllReduceLocality::SubscribeAllReduce =
        0;
int AllreduceAlgorithm::AllReducesInfo::AllReduceLocality::
    SubscribeAllReduce_handler(const hpx_t* allreduces, const size_t size) {
  NEUROX_MEM_PIN(uint64_t);
  assert(input_params_->allreduce_at_locality_);
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

hpx_action_t AllreduceAlgorithm::AllReducesInfo::AllReduceLocality::
    UnsubscribeAllReduce = 0;
int AllreduceAlgorithm::AllReducesInfo::AllReduceLocality::
    UnsubscribeAllReduce_handler(const hpx_t* allreduces, const size_t size) {
  NEUROX_MEM_PIN(uint64_t);
  assert(input_params_->allreduce_at_locality_);
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

hpx_action_t AllreduceAlgorithm::AllReducesInfo::Init = 0;
void AllreduceAlgorithm::AllReducesInfo::Init_handler(void*, const size_t) {}

hpx_action_t AllreduceAlgorithm::AllReducesInfo::Reduce = 0;
void AllreduceAlgorithm::AllReducesInfo::Reduce_handler(void*, const void*,
                                                        const size_t) {}

void AllreduceAlgorithm::AllReducesInfo::RegisterHpxActions() {
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
