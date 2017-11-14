#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::synchronizers;
using namespace neurox::interpolators;

Synchronizer* Synchronizer::New(SynchronizerIds type) {
  switch (type) {
    case SynchronizerIds::kDebug:
      return new DebugSynchronizer();
    case SynchronizerIds::kCoreneuron:
      return new CoreneuronSynchronizer();
    case SynchronizerIds::kAllReduce:
      return new AllreduceSynchronizer();
    case SynchronizerIds::kSlidingTimeWindow:
      return new SlidingTimeWindowSynchronizer();
    case SynchronizerIds::kTimeDependency:
      return new TimeDependencySynchronizer();
    default:
      return nullptr;
  }
  return nullptr;
};

SynchronizerNeuronInfo* SynchronizerNeuronInfo::New(SynchronizerIds type) {
  switch (type) {
    case SynchronizerIds::kDebug:
      return new DebugSynchronizer::CommunicationBarrier();
    case SynchronizerIds::kCoreneuron:
      return new CoreneuronSynchronizer::CommunicationBarrier();
    case SynchronizerIds::kAllReduce:
      return new AllreduceSynchronizer::AllReduceNeuronInfo(
          AllreduceSynchronizer::kAllReducesCount);
    case SynchronizerIds::kSlidingTimeWindow:
      // same as above, with different number of reduces
      return new AllreduceSynchronizer::AllReduceNeuronInfo(
          SlidingTimeWindowSynchronizer::kAllReducesCount);
    case SynchronizerIds::kTimeDependency:
      return new TimeDependencySynchronizer::TimeDependencies();
    default:
      return nullptr;
  }
  return nullptr;
}

hpx_action_t Synchronizer::CallInitLocality = 0;
int Synchronizer::CallInitLocality_handler(const int* synchronizer_id_ptr,
                                           const size_t) {
  NEUROX_MEM_PIN(uint64_t);

  // delete previous synchronizer (if any)
  delete neurox::synchronizer_;

  // initiate synchronizer
  SynchronizerIds synchronizer_id = *(SynchronizerIds*)synchronizer_id_ptr;
  neurox::synchronizer_ = Synchronizer::New(synchronizer_id);
  neurox::synchronizer_->InitLocality();

  // if we use "last neuron advances first" methodology
  if (input_params_->locality_comm_reduce_ &&
      synchronizer_->LocalitySyncInterval() == -1) {
    // scheduler semaphore (controls how many parallel jobs can run)
    size_t thread_count = hpx_get_num_threads();
    const int max_jobs = std::min(thread_count, locality::neurons_->size());
    locality::neurons_scheduler_sema_ = hpx_lco_sema_new(max_jobs);
    assert(max_jobs > 0);

    // neurons progress and its mutex
    locality::neurons_progress_mutex_ = hpx_lco_sema_new(1);
    locality::neurons_progress_ = new set<pair<floble_t, hpx_t>>();

    // create a trigger for each neuron, and broadcast it
    hpx_t lco = hpx_lco_and_new(neurox::locality::neurons_->size());
    const double tstart = input_params_->tstart_;
    for (hpx_t neuron_addr : *neurox::locality::neurons_) {
      // create an and-gate for stepping
      hpx_t step_trigger = hpx_lco_and_new(1);

      // inform neuron to add its information
      hpx_call(neuron_addr, Branch::SetSyncStepTrigger, lco, &step_trigger,
               sizeof(hpx_t));

      // set this value in the ordered set of neurons stepping
      locality::neurons_progress_->insert(std::make_pair(tstart, step_trigger));
    }
    hpx_lco_wait(lco);
    hpx_lco_delete(lco, HPX_NULL);
  }
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::NeuronInfoConstructor = 0;
int Synchronizer::NeuronInfoConstructor_handler(const int* synchronizer_id_ptr,
                                                const size_t) {
  NEUROX_MEM_PIN(Branch);
  SynchronizerIds sync_id = *(SynchronizerIds*)synchronizer_id_ptr;
  Neuron* soma = local->soma_;
  /* SynchronizerNeuronInfo may have been created before
   * by Neuron::Neuron(...) constructor */
  if (soma && soma->synchronizer_neuron_info_ == nullptr)
    local->soma_->synchronizer_neuron_info_ =
        SynchronizerNeuronInfo::New(sync_id);
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::CallInitNeuron = 0;
int Synchronizer::CallInitNeuron_handler() {
  NEUROX_MEM_PIN(Branch);
  NEUROX_RECURSIVE_BRANCH_ASYNC_CALL(Synchronizer::CallInitNeuron);
  assert(neurox::synchronizer_);
  neurox::synchronizer_->InitNeuron(local);
  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::RunLocality = 0;
int Synchronizer::RunLocality_handler(const double* tstop_ptr, const size_t) {
  NEUROX_MEM_PIN(uint64_t);

  if (locality::neurons_->empty()) NEUROX_MEM_UNPIN;

  const double tstop = *tstop_ptr;
  const double tstart = tstop - input_params_->tstop_;  // for all-benchmark

  const double reduction_dt = synchronizer_->LocalitySyncInterval();
  if (reduction_dt == 0)  // no reduction, launch neurons independently
  {
    wrappers::CallLocalNeurons(Synchronizer::RunNeuron, tstop_ptr,
                               sizeof(double));
  } else if (reduction_dt == -1)  // scheduler ON: step last-neuron first
  {
    // Launch neurons async. (they wait for the trigger to continue)
    const double tstop = *tstop_ptr;
    const size_t neurons_count = neurox::locality::neurons_->size();
    hpx_t neurons_lco = hpx_lco_and_new(neurons_count);

    for (size_t i = 0; i < neurons_count; i++)
      hpx_call(neurox::locality::neurons_->at(i), Synchronizer::RunNeuron,
               neurons_lco, tstop_ptr, sizeof(double));

    // number of simultaneous neuron to launch (1 per thread)
    auto neuron_it = locality::neurons_progress_->begin();
    double last_time = neuron_it->first;
    hpx_t trigger_addr = HPX_NULL;

    // TODO we must set all times of the pairs to 'tstart'!
    // otherwise when running benchmark it will start at 0,
    // instead of end time of previous synchronizer
    while (1) {
      // wait if all possilbe neurons are running
      hpx_lco_sema_p(locality::neurons_scheduler_sema_);

      // get trigger of last neuron and launch it
      //(delete it from set of jobs, will be added by neuron later)
      hpx_lco_sema_p(locality::neurons_progress_mutex_);
      neuron_it = locality::neurons_progress_->begin();
      last_time = neuron_it->first;
      trigger_addr = neuron_it->second;
      locality::neurons_progress_->erase(neuron_it);
      hpx_lco_sema_v_sync(locality::neurons_progress_mutex_);

      // if all jobs finished, exit
      if (last_time > tstop - 0.00001) break;

      // launch the last neuron in time
      hpx_lco_set(trigger_addr, 0, nullptr, HPX_NULL, HPX_NULL);
    }

    hpx_lco_wait_reset(neurons_lco);
    hpx_lco_delete_sync(neurons_lco);
  } else  // positive comm-reduce interval
  {
    double step_to_time = -1;
    for (double t = tstart; t <= tstop - 0.00001; t += reduction_dt) {
      synchronizer_->LocalitySyncInit();
      step_to_time = t + reduction_dt;
      wrappers::CallLocalNeurons(Synchronizer::RunNeuron, &step_to_time,
                                 sizeof(double));
      synchronizer_->LocalitySyncEnd();
    }
  }
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::RunNeuron = 0;
int Synchronizer::RunNeuron_handler(const double* tstop_ptr,
                                    const size_t size) {
  NEUROX_MEM_PIN(Branch);
  NEUROX_RECURSIVE_BRANCH_ASYNC_CALL(Synchronizer::RunNeuron, tstop_ptr, size);
  const double tstop = *tstop_ptr;
  Interpolator* interpolator = local->interpolator_;

  /* If not soma:
   * - launch job until tstop;
   * - Neuron-level reduction only needs to be performed by soma;
   * - Outgoing spikes handling only needs to be handled by soma;
   * - This branch will be synchronized with soma during Hines solver;
   */
  if (!local->soma_) {
    interpolator->StepTo(local, tstop);
    NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
    NEUROX_MEM_UNPIN;
  }

  /* If soma:
   * - it handles neuron-level reduction
   * - it waits for job scheduler (last-neuron first)
   * - it performs outgoing spikes handling
   */
  floble_t tpause = -1, max_step = -1;
  hpx_t spikes_lco = HPX_NULL;
  const hpx_t step_trigger = local->soma_->synchronizer_step_trigger_;
  const bool has_scheduler = step_trigger != HPX_NULL;
  double& t = local->nt_->_t;
  // const double dt_io = input_params_->dt_io_;

  while (t < tstop - 0.00001) {
    // do before-step operations e.g. mark step in all-reduces
    synchronizer_->NeuronSyncInit(local);

    // if scheduler is active: wait for scheduler signal to proceed
    if (has_scheduler) hpx_lco_wait_reset(step_trigger);

    // step to the next possible time instant or wait for one
    // if size too small or zero, wait to be awake again
    // (while it waits, allow scheduler to start a new job)
    max_step = synchronizer_->GetNeuronMaxStep(local);
    assert(max_step >= 0.025);

    tpause = std::min(t + max_step, tstop);
    assert(tpause > t + 0.000001);  // make sure it will step
    spikes_lco = interpolator->StepTo(local, tpause);
    assert(fabs(t - tpause) < 0.000001);  // equal

    // do after-step operations e.g. wait for spike delivery
    synchronizer_->NeuronSyncEnd(local, spikes_lco);
    // if (fmod(t, dt_io) == 0) {  /*output*/ }

    // decrement scheduler semaphore (wake up if necessary)
    if (has_scheduler) {
      // increment scheduler counter to allow it to look for next job
      hpx_lco_sema_v_sync(locality::neurons_scheduler_sema_);

      // re-add this job to queue, to be picked up again later
      hpx_lco_sema_p(locality::neurons_progress_mutex_);
      // hack: make it not be picked by scheduler immediately
      tpause += 0.00000001;
      locality::neurons_progress_->insert(std::make_pair(tpause, step_trigger));
      hpx_lco_sema_v_sync(locality::neurons_progress_mutex_);
    } else {
      // wait for time dependencies!
    }
  }
  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
#ifndef NDEBUG
  printf("## Neuron %d finished\n", local->soma_->gid_);
#endif
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::CallClearLocality = 0;
int Synchronizer::CallClearLocality_handler() {
  NEUROX_MEM_PIN(uint64_t);

  if (input_params_->locality_comm_reduce_ &&
      synchronizer_->LocalitySyncInterval() == -1) {
    hpx_lco_delete_sync(neurox::locality::neurons_progress_mutex_);
    hpx_lco_delete_sync(neurox::locality::neurons_scheduler_sema_);

    for (auto& neuron_it : *locality::neurons_progress_)
      hpx_lco_delete_sync(neuron_it.second);

    // set step trigger to HPX_NULL in all branches
    hpx_t step_trigger = HPX_NULL;
    hpx_t lco = hpx_lco_and_new(neurox::locality::neurons_->size());
    for (hpx_t neuron_addr : *neurox::locality::neurons_)
      hpx_call(neuron_addr, Branch::SetSyncStepTrigger, lco, &step_trigger,
               sizeof(hpx_t));
    hpx_lco_wait(lco);
    hpx_lco_delete(lco, HPX_NULL);

    (*neurox::locality::neurons_progress_).clear();
    delete neurox::locality::neurons_progress_;
    neurox::locality::neurons_progress_ = nullptr;
  }
  neurox::synchronizer_->ClearLocality();
  delete neurox::synchronizer_;
  neurox::synchronizer_ = nullptr;
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::CallClearNeuron = 0;
int Synchronizer::CallClearNeuron_handler() {
  NEUROX_MEM_PIN(Branch);
  NEUROX_RECURSIVE_BRANCH_ASYNC_CALL(Synchronizer::CallClearNeuron);
  neurox::synchronizer_->ClearNeuron(local);
  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::NeuronInfoDestructor = 0;
int Synchronizer::NeuronInfoDestructor_handler() {
  NEUROX_MEM_PIN(Branch);
  if (local->soma_) {
    delete local->soma_->synchronizer_neuron_info_;
    local->soma_->synchronizer_neuron_info_ = nullptr;
  }
  NEUROX_MEM_UNPIN;
}

/// auxiliar method for CallLocalNeurons
hpx_action_t Synchronizer::CallAllNeuronsAux = 0;
int Synchronizer::CallAllNeuronsAux_handler(const int nargs, const void* args[],
                                            const size_t sizes[]) {
  // TODO clean this
  NEUROX_MEM_PIN(uint64_t);
  hpx_action_t f = *(hpx_action_t*)args[0];  // first arg is action id
  if (nargs == 1)
    wrappers::CallLocalNeurons(f);
  else if (nargs == 2)
    wrappers::CallLocalNeurons(f, args[1], sizes[1]);
  else if (nargs == 3)
    wrappers::CallLocalNeurons(f, args[1], sizes[1], args[2], sizes[2]);
  else if (nargs == 4)
    wrappers::CallLocalNeurons(f, args[1], sizes[1], args[2], sizes[2], args[3],
                               sizes[3]);
  else if (nargs == 5)
    wrappers::CallLocalNeurons(f, args[1], sizes[1], args[2], sizes[2], args[3],
                               sizes[3], args[4], sizes[4]);
  else if (nargs == 6)
    wrappers::CallLocalNeurons(f, args[1], sizes[1], args[2], sizes[2], args[3],
                               sizes[3], args[4], sizes[4], args[5], sizes[5]);
  else {
    assert(0);
  }
  NEUROX_MEM_UNPIN;
}

void Synchronizer::RegisterHpxActions() {
  wrappers::RegisterMultipleVarAction(Synchronizer::CallAllNeuronsAux,
                                      Synchronizer::CallAllNeuronsAux_handler);

  wrappers::RegisterZeroVarAction(Synchronizer::CallClearLocality,
                                  Synchronizer::CallClearLocality_handler);
  wrappers::RegisterZeroVarAction(Synchronizer::CallClearNeuron,
                                  Synchronizer::CallClearNeuron_handler);
  wrappers::RegisterZeroVarAction(Synchronizer::CallInitNeuron,
                                  Synchronizer::CallInitNeuron_handler);
  wrappers::RegisterZeroVarAction(Synchronizer::NeuronInfoDestructor,
                                  Synchronizer::NeuronInfoDestructor_handler);
  wrappers::RegisterSingleVarAction<double>(Synchronizer::RunLocality,
                                            Synchronizer::RunLocality_handler);
  wrappers::RegisterSingleVarAction<double>(Synchronizer::RunNeuron,
                                            Synchronizer::RunNeuron_handler);
  wrappers::RegisterSingleVarAction<int>(
      Synchronizer::CallInitLocality, Synchronizer::CallInitLocality_handler);
  wrappers::RegisterSingleVarAction<int>(
      Synchronizer::NeuronInfoConstructor,
      Synchronizer::NeuronInfoConstructor_handler);

  // registration of instantiated synchronizers
  AllreduceSynchronizer::RegisterHpxActions();
  TimeDependencySynchronizer::RegisterHpxActions();
}
