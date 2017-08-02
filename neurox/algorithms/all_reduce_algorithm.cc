#include "neurox/algorithms/all_reduce_algorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

hpx_t* AllReduceAlgorithm::allReduces = nullptr;

AllReduceAlgorithm::AllReduceAlgorithm() {
  AllReduceAlgorithm::AllReducesInfo::reductionsPerCommStep =
      AllReduceAlgorithm::kAllReducesCount;
}

AllReduceAlgorithm::~AllReduceAlgorithm() {}

const AlgorithmType AllReduceAlgorithm::GetType() {
  return AlgorithmType::kBackwardEulerAllReduce;
}

const char* AllReduceAlgorithm::GetTypeString() {
  return "BackwardEulerAllReduce";
}

void AllReduceAlgorithm::Init() {
  SubscribeAllReduces(AllReduceAlgorithm::allReduces,
                      AllReduceAlgorithm::kAllReducesCount);
}

void AllReduceAlgorithm::Clear() {
  UnsubscribeAllReduces(AllReduceAlgorithm::allReduces,
                        AllReduceAlgorithm::kAllReducesCount);
}

double AllReduceAlgorithm::Launch() {
  int totalSteps = Algorithm::getTotalStepsCount();
  hpx_time_t now = hpx_time_now();
  if (input_params->all_reduce_at_locality_)
    hpx_bcast_rsync(Branch::BackwardEulerOnLocality, &totalSteps, sizeof(int));
  else
    neurox::wrappers::CallAllNeurons(Branch::BackwardEuler, &totalSteps,
                                     sizeof(int));
  double elapsedTime = hpx_time_elapsed_ms(now) / 1e3;
  input::Debugger::RunCoreneuronAndCompareAllBranches();
  return elapsedTime;
}

void AllReduceAlgorithm::StepBegin(Branch*) {}

void AllReduceAlgorithm::StepEnd(Branch* b, hpx_t spikesLco) {
  WaitForSpikesDelivery(b, spikesLco);
  input::Debugger::SingleNeuronStepAndCompare(&nrn_threads[b->nt_->id], b,
                                              input_params->second_order_);
}

void AllReduceAlgorithm::Run(Branch* b, const void* args) { Run2(b, args); }

hpx_t AllReduceAlgorithm::SendSpikes(Neuron* b, double tt, double) {
  return SendSpikes2(b, tt);
}

void AllReduceAlgorithm::SubscribeAllReduces(hpx_t*& allReduces,
                                             size_t allReducesCount) {
  assert(allReduces == nullptr);
  allReduces = new hpx_t[allReducesCount];

  hpx_bcast_rsync(AllReduceAlgorithm::AllReducesInfo::SetReductionsPerCommStep,
                  &allReducesCount, sizeof(int));

  for (int i = 0; i < allReducesCount; i++)
    allReduces[i] = hpx_process_collective_allreduce_new(
        0, AllReduceAlgorithm::AllReducesInfo::Init,
        AllReduceAlgorithm::AllReducesInfo::Reduce);

  if (input_params->all_reduce_at_locality_)
    hpx_bcast_rsync(AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::
                        SubscribeAllReduce,
                    allReduces, sizeof(hpx_t) * allReducesCount);
  else
    neurox::wrappers::CallAllNeurons(
        AllReduceAlgorithm::AllReducesInfo::SubscribeAllReduce, allReduces,
        sizeof(hpx_t) * allReducesCount);

  for (int i = 0; i < allReducesCount; i++)
    hpx_process_collective_allreduce_subscribe_finalize(allReduces[i]);
}

void AllReduceAlgorithm::UnsubscribeAllReduces(hpx_t*& allReduces,
                                               size_t allReducesCount) {
  assert(allReduces != nullptr);
  if (input_params->all_reduce_at_locality_)
    hpx_bcast_rsync(AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::
                        UnsubscribeAllReduce,
                    allReduces, sizeof(hpx_t) * allReducesCount);
  else
    neurox::wrappers::CallAllNeurons(
        AllReduceAlgorithm::AllReducesInfo::UnsubscribeAllReduce, allReduces,
        sizeof(hpx_t) * allReducesCount);

  for (int i = 0; i < allReducesCount; i++)
    hpx_process_collective_allreduce_delete(allReduces[i]);

  delete[] allReduces;
  allReduces = nullptr;
}

void AllReduceAlgorithm::WaitForSpikesDelivery(Branch* b, hpx_t spikesLco) {
  // wait for spikes sent 4 steps ago (queue has always size 3)
  if (b->soma_) {
    AllReducesInfo* stw = (AllReducesInfo*)b->soma_->algorithm_metadata_;
    std::queue<hpx_t> q = stw->spikesLcoQueue;
    assert(stw->spikesLcoQueue.size() ==
           CoreneuronAlgorithm::CommunicationBarrier::kCommStepSize - 1);
    stw->spikesLcoQueue.push(spikesLco);
    hpx_t queuedSpikesLco = stw->spikesLcoQueue.front();
    stw->spikesLcoQueue.pop();
    if (queuedSpikesLco != HPX_NULL) {
      hpx_lco_wait(queuedSpikesLco);
      hpx_lco_delete_sync(queuedSpikesLco);
    }
  }
}

hpx_t AllReduceAlgorithm::SendSpikes2(Neuron* neuron, double tt) {
  hpx_t newSynapsesLco = hpx_lco_and_new(neuron->synapses_.size());
  for (Neuron::Synapse*& s : neuron->synapses_)
    hpx_call(s->branch_addr_, Branch::AddSpikeEvent, newSynapsesLco, &neuron->gid_,
             sizeof(neuron_id_t), &tt, sizeof(spike_time_t));
  return newSynapsesLco;
}

void AllReduceAlgorithm::Run2(Branch* b, const void* args) {
  int steps = *(int*)args;
  const int reductionsPerCommStep =
      AllReduceAlgorithm::AllReducesInfo::reductionsPerCommStep;
  const int commStepSize =
      CoreneuronAlgorithm::CommunicationBarrier::kCommStepSize;
  const int stepsPerReduction = commStepSize / reductionsPerCommStep;
  const AllReducesInfo* stw =
      b->soma_ ? (AllReducesInfo*)b->soma_->algorithm_metadata_ : nullptr;

  for (int s = 0; s < steps; s += commStepSize)  // for every communication step
  {
#ifndef NDEBUG
    if (hpx_get_my_rank() == 0 && b->nt_->id == 0 && b->soma_) {
      printf("-- t=%.4f ms\n", input_params->dt_ * s);
      fflush(stdout);
    }
#endif
    for (int r = 0; r < reductionsPerCommStep; r++)  // for every reduction step
    {
      if (b->soma_) {
        if (s >= commStepSize)  // first comm-window does not wait
          hpx_lco_wait_reset(stw->allReduceFuture[r]);
        else
          // fixes crash for Algorithm::ALL when running two hpx-reduce -based
          // algorithms in a row
          hpx_lco_reset_sync(stw->allReduceFuture[r]);

        hpx_process_collective_allreduce_join(stw->allReduceLco[r],
                                              stw->allReduceId[r], NULL, 0);
      }

      for (int n = 0; n < stepsPerReduction; n++) b->BackwardEulerStep();
      // Input::Coreneuron::Debugger::stepAfterStepBackwardEuler(local,
      // &nrn_threads[this->nt->id], secondorder); //SMP ONLY
    }
  }
}

AllReduceAlgorithm::AllReducesInfo::AllReducesInfo() {
  for (int s = 0;
       s < CoreneuronAlgorithm::CommunicationBarrier::kCommStepSize - 1; s++)
    this->spikesLcoQueue.push(HPX_NULL);
}

AllReduceAlgorithm::AllReducesInfo::~AllReducesInfo() {
  for (int i = 0; i < spikesLcoQueue.size(); i++) {
    hpx_t queuedSpikesLco = spikesLcoQueue.front();
    if (queuedSpikesLco != HPX_NULL) hpx_lco_delete_sync(queuedSpikesLco);
    spikesLcoQueue.pop();
  }
}

hpx_action_t AllReduceAlgorithm::AllReducesInfo::SubscribeAllReduce = 0;
int AllReduceAlgorithm::AllReducesInfo::SubscribeAllReduce_handler(
    const hpx_t* allreduces, const size_t size) {
  NEUROX_MEM_PIN(Branch);
  AllReducesInfo* stw = (AllReducesInfo*)local->soma_->algorithm_metadata_;
  stw->allReduceFuture = new hpx_t[AllReducesInfo::reductionsPerCommStep];
  stw->allReduceLco = new hpx_t[AllReducesInfo::reductionsPerCommStep];
  stw->allReduceId = new int[AllReducesInfo::reductionsPerCommStep];
  for (int i = 0; i < size / sizeof(hpx_t); i++) {
    stw->allReduceLco[i] = allreduces[i];
    stw->allReduceFuture[i] = hpx_lco_future_new(0);  // no value to be reduced
    stw->allReduceId[i] = hpx_process_collective_allreduce_subscribe(
        allreduces[i], hpx_lco_set_action, stw->allReduceFuture[i]);
  }
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t AllReduceAlgorithm::AllReducesInfo::UnsubscribeAllReduce = 0;
int AllReduceAlgorithm::AllReducesInfo::UnsubscribeAllReduce_handler(
    const hpx_t* allreduces, const size_t size) {
  NEUROX_MEM_PIN(Branch);
  AllReducesInfo* stw = (AllReducesInfo*)local->soma_->algorithm_metadata_;
  for (int i = 0; i < size / sizeof(hpx_t); i++) {
    hpx_process_collective_allreduce_unsubscribe(allreduces[i],
                                                 stw->allReduceId[i]);
    if (stw->allReduceFuture[i] != HPX_NULL)
      hpx_lco_delete_sync(stw->allReduceFuture[i]);
  }
  delete[] stw->allReduceLco;
  stw->allReduceLco = nullptr;
  delete[] stw->allReduceFuture;
  stw->allReduceFuture = nullptr;
  delete[] stw->allReduceId;
  stw->allReduceId = nullptr;
  return neurox::wrappers::MemoryUnpin(target);
}

int AllReduceAlgorithm::AllReducesInfo::reductionsPerCommStep = -1;
std::vector<hpx_t>*
    AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::localityNeurons =
        nullptr;
hpx_t* AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::allReduceFuture =
    nullptr;
hpx_t* AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::allReduceLco =
    nullptr;
int* AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::allReduceId =
    nullptr;

hpx_action_t AllReduceAlgorithm::AllReducesInfo::SetReductionsPerCommStep = 0;
int AllReduceAlgorithm::AllReducesInfo::SetReductionsPerCommStep_handler(
    const int* val, const size_t) {
  NEUROX_MEM_PIN(uint64_t);
  reductionsPerCommStep = *val;
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t
    AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::SubscribeAllReduce =
        0;
int AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::
    SubscribeAllReduce_handler(const hpx_t* allreduces, const size_t size) {
  NEUROX_MEM_PIN(uint64_t);
  assert(input_params->all_reduce_at_locality_);
  AllReduceLocality::allReduceLco =
      new hpx_t[AllReducesInfo::reductionsPerCommStep];
  AllReduceLocality::allReduceFuture =
      new hpx_t[AllReducesInfo::reductionsPerCommStep];
  AllReduceLocality::allReduceId =
      new int[AllReducesInfo::reductionsPerCommStep];
  for (int i = 0; i < size / sizeof(hpx_t); i++) {
    allReduceLco[i] = allreduces[i];
    allReduceFuture[i] = hpx_lco_future_new(0);  // no value to be reduced
    allReduceId[i] = hpx_process_collective_allreduce_subscribe(
        allreduces[i], hpx_lco_set_action, allReduceFuture[i]);
  }
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::
    UnsubscribeAllReduce = 0;
int AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::
    UnsubscribeAllReduce_handler(const hpx_t* allreduces, const size_t size) {
  NEUROX_MEM_PIN(uint64_t);
  assert(input_params->all_reduce_at_locality_);
  for (int i = 0; i < size / sizeof(hpx_t); i++) {
    hpx_process_collective_allreduce_unsubscribe(allreduces[i], allReduceId[i]);
    hpx_lco_delete_sync(allReduceFuture[i]);
  }
  delete[] allReduceLco;
  allReduceLco = nullptr;
  delete[] allReduceFuture;
  allReduceFuture = nullptr;
  delete[] allReduceId;
  allReduceId = nullptr;
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t AllReduceAlgorithm::AllReducesInfo::Init = 0;
void AllReduceAlgorithm::AllReducesInfo::Init_handler(void*, const size_t) {}

hpx_action_t AllReduceAlgorithm::AllReducesInfo::Reduce = 0;
void AllReduceAlgorithm::AllReducesInfo::Reduce_handler(void*, const void*,
                                                        const size_t) {}

void AllReduceAlgorithm::AllReducesInfo::RegisterHpxActions() {
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
