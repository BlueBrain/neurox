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

  if (input_params_->output_comm_count_)
    Statistics::CommCount::mutex = hpx_lco_sema_new(1);

  // delete previous synchronizer (if any)
  delete neurox::synchronizer_;

  // initiate synchronizer
  SynchronizerIds synchronizer_id = *(SynchronizerIds*)synchronizer_id_ptr;
  neurox::synchronizer_ = Synchronizer::New(synchronizer_id);
  neurox::synchronizer_->InitLocality();

  // if we use "last neuron advances first" methodology
  if (  // input_params_->locality_comm_reduce_ ||
      input_params_->scheduler_) {
    // scheduler semaphore (controls how many parallel jobs can run)
    size_t thread_count = hpx_get_num_threads();
    const int max_jobs = std::min(thread_count, locality::neurons_->size());
    assert(max_jobs > 0);
    locality::scheduler_neurons_sema_ = hpx_lco_sema_new(max_jobs);
#if defined(PRINT_TIME_DEPENDENCY) or defined(PRINT_TIME_DEPENDENCY_MUTEX) or \
    defined(PRINT_TIME_DEPENDENCY_STEP_SIZE)
    locality::scheduler_sema_counter_ = max_jobs;
#endif

    // neurons progress and its mutex
    libhpx_cond_init(&locality::scheduler_wait_condition_);
    libhpx_mutex_init(&locality::scheduler_lock_);

    locality::scheduler_neurons_ = new set<pair<floble_t, hpx_t>>();

    // create a trigger for each neuron, and broadcast it
    hpx_t lco = hpx_lco_and_new(neurox::locality::neurons_->size());
    const double tstart = input_params_->tstart_;
    for (hpx_t neuron_addr : *neurox::locality::neurons_) {
      hpx_t step_trigger = hpx_lco_and_new(1);
      hpx_call(neuron_addr, Branch::SetSyncStepTrigger, lco, &step_trigger,
               sizeof(hpx_t));

      // set this value in the set of next neurons stepping
      locality::scheduler_neurons_->insert(make_pair(tstart, step_trigger));
    }
    hpx_lco_wait(lco);
    hpx_lco_delete(lco, HPX_NULL);
  }
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

  if (locality::neurons_->empty()) NEUROX_MEM_UNPIN;  // no neurons... no run!

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
    locality::scheduler_remaining_neurons_ = locality::neurons_->size();
    hpx_t neurons_lco = hpx_lco_and_new(locality::scheduler_remaining_neurons_);
    for (size_t i = 0; i < locality::scheduler_remaining_neurons_; i++)
      hpx_call(neurox::locality::neurons_->at(i), Synchronizer::RunNeuron,
               neurons_lco, tstop_ptr, sizeof(double));

    libhpx_mutex_lock(&locality::scheduler_lock_);
    while (locality::scheduler_remaining_neurons_ > 1) {
      libhpx_mutex_unlock(&locality::scheduler_lock_);

#ifdef PRINT_TIME_DEPENDENCY_MUTEX
      fprintf(stderr,
              "## %d.%d.%d ## before "
              "hpx_lco_sema(locality::neurons_scheduler_sema_ %d) 1\n",
              wrappers::MyRankId(), wrappers::MyThreadId(),
              wrappers::MyLightWeightThreadId(),
              locality::scheduler_sema_counter_);
#endif
      hpx_lco_sema_p(locality::scheduler_neurons_sema_);
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
      fprintf(stderr,
              "## %d.%d.%d ## after  "
              "hpx_lco_sema(locality::neurons_scheduler_sema_ %d) 1\n",
              wrappers::MyRankId(), wrappers::MyThreadId(),
              wrappers::MyLightWeightThreadId(),
              --locality::scheduler_sema_counter_);
#endif

#ifdef PRINT_TIME_DEPENDENCY_MUTEX
      fprintf(stderr,
              "## %d.%d ## before "
              "libhpx_mutex_lock(&locality::scheduler_lock_) 2\n",
              wrappers::MyRankId(), wrappers::MyThreadId());
#endif
      // mark this thread sleep until there's an entry in set
      // or until one neuron has finished (search libhpx_cond_broadcast)
      libhpx_mutex_lock(&locality::scheduler_lock_);
      if (locality::scheduler_neurons_->empty())
        libhpx_cond_wait(&locality::scheduler_wait_condition_,
                         &locality::scheduler_lock_);
      if (!locality::scheduler_neurons_->empty()) {
        auto next_neuron = locality::scheduler_neurons_->begin();
        locality::scheduler_neurons_->erase(next_neuron);

        // launch the last neuron in time or queued
        hpx_lco_set(next_neuron->second, 0, nullptr, HPX_NULL, HPX_NULL);
#ifdef PRINT_TIME_DEPENDENCY
        fprintf(stderr, "## %d.%d ## scheduler adds to queue %d %.12f\n",
                wrappers::MyRankId(), wrappers::MyThreadId(),
                locality::from_hpx_to_gid->at(next_neuron->second),
                next_neuron->first);
#endif
      }
    }
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
    fprintf(stderr, "## %d.%d.%d ## scheduler waiting for neurons_lco 99\n",
            wrappers::MyRankId(), wrappers::MyThreadId(),
            wrappers::MyLightWeightThreadId());
#endif
    libhpx_mutex_unlock(&locality::scheduler_lock_);
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
    BackwardEuler::Step(local);
    interpolator->StepTo(local, tstop);
    NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
    NEUROX_MEM_UNPIN;
  }

  /* If soma:
   * - it handles neuron-level reduction
   * - it waits for job scheduler (last-neuron first)
   * - it performs outgoing spikes handling
   */
  floble_t tpause = -1, dt_pause = -1;
  const hpx_t step_trigger = local->soma_->scheduler_step_trigger_;
  const bool has_scheduler = step_trigger != HPX_NULL;
  NrnThread* nt = local->nt_;
  double& t = nt->_t;
  // const double dt_io = input_params_->dt_io_;

  while (t < tstop - 1e-5) {
    // do before-step operations e.g. mark step in all-reduces
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
    fprintf(stderr, "~~ %d.%d ~~ before NeuronSyncInit 3\n",
            wrappers::MyRankId(), wrappers::MyThreadId());
#endif
    synchronizer_->NeuronSyncInit(local);

#ifdef PRINT_TIME_DEPENDENCY_MUTEX
    fprintf(stderr,
            "~~ %d.%d.%d ~~ before "
            "hpx_lco_sema(locality::neurons_scheduler_sema_ %d) 4\n",
            wrappers::MyRankId(), wrappers::MyThreadId(),
            wrappers::MyLightWeightThreadId(),
            locality::scheduler_sema_counter_);
#endif
    // if scheduler is active: wait for scheduler signal to proceed
    if (has_scheduler) hpx_lco_wait_reset(step_trigger);

    // step to the next possible time instant or wait for one if not scheduler
    dt_pause = synchronizer_->GetNeuronMaxStep(local);
    tpause = std::min(t + dt_pause, tstop) + 1e-8;

#ifdef PRINT_TIME_DEPENDENCY_STEP_SIZE
    if (has_scheduler)
      fprintf(stderr, "step,%d,%d,%.4f,%.4f,%.4f\n", neurox::neurons_count_,
              local->soma_->gid_, t, tpause, tpause - t);
#endif
    interpolator->StepTo(local, tpause);

    // do after-step operations e.g. wait for spike delivery
    synchronizer_->NeuronSyncEnd(local);
    // if (fmod(t, dt_io) == 0) {  /*output*/ }
  }
  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::CallClearLocality = 0;
int Synchronizer::CallClearLocality_handler() {
  NEUROX_MEM_PIN(uint64_t);

  if (input_params_->output_comm_count_)
    hpx_lco_delete_sync(Statistics::CommCount::mutex);

  if (input_params_->locality_comm_reduce_ &&
      synchronizer_->LocalitySyncInterval() == -1) {
    libhpx_cond_destroy(&locality::scheduler_wait_condition_);
    libhpx_mutex_destroy(&locality::scheduler_lock_);
    hpx_lco_delete_sync(locality::scheduler_neurons_sema_);

    for (auto& neuron_it : *locality::scheduler_neurons_)
      hpx_lco_delete_sync(neuron_it.second);

    // set step trigger to HPX_NULL in all branches
    hpx_t step_trigger = HPX_NULL;
    hpx_t lco = hpx_lco_and_new(neurox::locality::neurons_->size());
    for (hpx_t neuron_addr : *neurox::locality::neurons_)
      hpx_call(neuron_addr, Branch::SetSyncStepTrigger, lco, &step_trigger,
               sizeof(hpx_t));
    hpx_lco_wait(lco);
    hpx_lco_delete(lco, HPX_NULL);

    (*neurox::locality::scheduler_neurons_).clear();
    delete neurox::locality::scheduler_neurons_;
    neurox::locality::scheduler_neurons_ = nullptr;
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
  wrappers::RegisterSingleVarAction<double>(Synchronizer::RunLocality,
                                            Synchronizer::RunLocality_handler);
  wrappers::RegisterSingleVarAction<double>(Synchronizer::RunNeuron,
                                            Synchronizer::RunNeuron_handler);
  wrappers::RegisterSingleVarAction<int>(
      Synchronizer::CallInitLocality, Synchronizer::CallInitLocality_handler);

  // registration of instantiated synchronizers
  AllreduceSynchronizer::RegisterHpxActions();
  TimeDependencySynchronizer::RegisterHpxActions();
}
