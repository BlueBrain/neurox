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
  if (input_params_->locality_comm_reduce_ ||
      input_params_->neurons_scheduler_) {
    // scheduler semaphore (controls how many parallel jobs can run)
    size_t thread_count = hpx_get_num_threads();
    // Note: i tried thread_count*2, not good: it does not process
    // pending messages RPCs
    const int max_jobs = std::min(thread_count, locality::neurons_->size());
    locality::neurons_scheduler_sema_ = hpx_lco_sema_new(max_jobs);
    assert(max_jobs > 0);

    // neurons progress and its mutex
    locality::neurons_progress_mutex_ = hpx_lco_sema_new(1);
    locality::neurons_progress_ = new set<pair<floble_t, hpx_t>>();
    locality::neurons_progress_queue_ = new std::queue<hpx_t>();

    // create a trigger for each neuron, and broadcast it
    hpx_t lco = hpx_lco_and_new(neurox::locality::neurons_->size());
    const double tstart = input_params_->tstart_;
    for (hpx_t neuron_addr : *neurox::locality::neurons_) {
      // create an and-gate for stepping
      hpx_t step_trigger = hpx_lco_and_new(1);

      // inform neuron to add its information
      hpx_call(neuron_addr, Branch::SetSyncStepTrigger, lco, &step_trigger,
               sizeof(hpx_t));

      // set this value in the queue of next neurons stepping
      locality::neurons_progress_queue_->push(step_trigger);

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

  if (locality::neurons_->empty()) NEUROX_MEM_UNPIN; //no neurons... no run!

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
    unsigned int remaining_neurons_count = locality::neurons_->size();
    hpx_t neurons_lco = hpx_lco_and_new(remaining_neurons_count);
    for (size_t i = 0; i < remaining_neurons_count; i++)
      hpx_call(neurox::locality::neurons_->at(i), Synchronizer::RunNeuron,
               neurons_lco, tstop_ptr, sizeof(double));

    // While last neuron is not done (ie all done)
    while (remaining_neurons_count>0) {
      // wait if all allowed neurons are running
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
      fprintf(stderr, "## %d.%d.%d ## before hpx_lco_sema(locality::neurons_scheduler_sema_) 1\n", 
              wrappers::MyRankId(), wrappers::MyThreadId(), wrappers::MyLightWeightThreadId());
#endif
      hpx_lco_sema_p(locality::neurons_scheduler_sema_);
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
      fprintf(stderr, "## %d.%d.%d ## after  hpx_lco_sema(locality::neurons_scheduler_sema_) 1\n", 
              wrappers::MyRankId(), wrappers::MyThreadId(), wrappers::MyLightWeightThreadId());
#endif

/*
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
      fprintf(stderr, "## %d.%d ## before hpx_lco_sema(locality::neurons_progress_mutex_) 2\n", 
              wrappers::MyRankId(), wrappers::MyThreadId());
#endif
*/
      hpx_lco_sema_p(locality::neurons_progress_mutex_);
/*
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
      fprintf(stderr, "## %d.%d ## after  hpx_lco_sema(locality::neurons_progress_mutex_) 2\n", 
              wrappers::MyRankId(), wrappers::MyThreadId());
#endif
*/
      //if there are no neurons queued to be run next, get them from the set of time ordered neurons
      if (locality::neurons_progress_queue_->empty()) {

        //add to queue all neurons with similar time instant
        auto last_neuron = locality::neurons_progress_->begin();
        while (!locality::neurons_progress_->empty() && std::fabs(locality::neurons_progress_->begin()->first - last_neuron->first) < 0.025) {
            //TODO increase this 0.25 value?!
            last_neuron = locality::neurons_progress_->begin();
            locality::neurons_progress_->erase(last_neuron);

            if (last_neuron->first >= tstop)
            {
              remaining_neurons_count--;
#ifdef PRINT_TIME_DEPENDENCY
              fprintf(stderr, "## %d.%d.%d ## Neuron %d finished. Remaining %d.\n", 
                      wrappers::MyRankId(), wrappers::MyThreadId(), wrappers::MyLightWeightThreadId(),
                      locality::from_hpx_to_gid->at(last_neuron->second), remaining_neurons_count);
#endif
            }
            else
            {
              locality::neurons_progress_queue_->push(last_neuron->second);
#ifdef PRINT_TIME_DEPENDENCY
              fprintf(stderr, "## %d.%d ## scheduler adds to queue %d %.12f\n",
                      wrappers::MyRankId(), wrappers::MyThreadId(),
                      locality::from_hpx_to_gid->at(last_neuron->second),
                      last_neuron->first);
#endif
            }
          }
      }

      if (!locality::neurons_progress_queue_->empty())
      {
        hpx_t last_neuron_addr = locality::neurons_progress_queue_->front();
        locality::neurons_progress_queue_->pop();
#ifdef PRINT_TIME_DEPENDENCY
        fprintf(stderr, "## %d.%d ## scheduler launches from queue %d (size queue %d)\n",
                wrappers::MyRankId(), wrappers::MyThreadId(),
                locality::from_hpx_to_gid->at(last_neuron_addr),
              locality::neurons_progress_queue_->size());
#endif
        // launch the last neuron in time or queued
        hpx_lco_set(last_neuron_addr, 0, nullptr, HPX_NULL, HPX_NULL);
      }
/*
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
      fprintf(stderr, "## %d.%d ## before hpx_lco_sema(locality::neurons_progress_mutex_) 3\n", 
              wrappers::MyRankId(), wrappers::MyThreadId());
#endif
*/
      hpx_lco_sema_v_sync(locality::neurons_progress_mutex_);
/*
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
      fprintf(stderr, "## %d.%d ## after hpx_lco_sema(locality::neurons_progress_mutex_) 3\n", 
              wrappers::MyRankId(), wrappers::MyThreadId());
#endif
*/
    }
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
    fprintf(stderr, "## %d.%d.%d ## scheduler waiting for neurons_lco 99\n", 
              wrappers::MyRankId(), wrappers::MyThreadId(), wrappers::MyLightWeightThreadId());
#endif
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
  floble_t tpause = -1, dtp = -1;
  hpx_t spikes_lco = HPX_NULL;
  const hpx_t step_trigger = local->soma_->synchronizer_step_trigger_;
  const bool has_scheduler = step_trigger != HPX_NULL;
  NrnThread* nt = local->nt_;
  double& t = nt->_t;
  // const double dt_io = input_params_->dt_io_;

  while (t < tstop - 0.00001) {
    // do before-step operations e.g. mark step in all-reduces
    synchronizer_->NeuronSyncInit(local);

    // if scheduler is active: wait for scheduler signal to proceed
    if (has_scheduler) hpx_lco_wait_reset(step_trigger);

    // step to the next possible time instant or wait for one
    // if size too small or zero, wait to be awake again
    // (while it waits, allow scheduler to start a new job)
    dtp = synchronizer_->GetNeuronMaxStep(local);
    tpause = std::min(t + dtp, tstop);

    // It may happen that notifications from dependencies haven't
    // arrived yet. In that case, go back to the queue of next neurons
    // to be picked up
    if (has_scheduler && tpause < t + 0.000001)  // if step too small or zero
    {
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
      fprintf(stderr, "~~ %d.%d.%d ~~ before hpx_lco_sema(locality::neurons_scheduler_sema_) 4\n", 
              wrappers::MyRankId(), wrappers::MyThreadId(), wrappers::MyLightWeightThreadId());
#endif
      hpx_lco_sema_v_sync(locality::neurons_scheduler_sema_);
      //hpx_lco_sema_v(locality::neurons_scheduler_sema_, HPX_NULL);
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
      fprintf(stderr, "~~ %d.%d.%d ~~ after hpx_lco_sema(locality::neurons_scheduler_sema_) 4\n", 
              wrappers::MyRankId(), wrappers::MyThreadId(),wrappers::MyLightWeightThreadId());
#endif
/*
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
      fprintf(stderr, "~~ %d.%d ~~ before hpx_lco_sema(locality::neurons_progress_mutex_) 5\n", 
              wrappers::MyRankId(), wrappers::MyThreadId());
#endif
*/
      hpx_lco_sema_p(locality::neurons_progress_mutex_);
/*
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
      fprintf(stderr, "~~ %d.%d ~~ after hpx_lco_sema(locality::neurons_progress_mutex_) 5\n", 
              wrappers::MyRankId(), wrappers::MyThreadId());
#endif
*/
      locality::neurons_progress_queue_->push(step_trigger);
/*
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
      fprintf(stderr, "~~ %d.%d ~~ before hpx_lco_sema(locality::neurons_progress_mutex_) 6\n", 
              wrappers::MyRankId(), wrappers::MyThreadId());
#endif
*/
      hpx_lco_sema_v_sync(locality::neurons_progress_mutex_);
/*
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
      fprintf(stderr, "~~ %d.%d ~~ after hpx_lco_sema(locality::neurons_progress_mutex_) 6\n", 
              wrappers::MyRankId(), wrappers::MyThreadId());
#endif
*/
    } else {
#ifdef PRINT_TIME_DEPENDENCY_STEP_SIZE
      if (has_scheduler)
        fprintf(stderr, "step,%d,%d,%.4f,%.4f,%.4f\n", neurox::neurons_count_,
                local->soma_->gid_, t, tpause, tpause - t);
#endif
      spikes_lco = interpolator->StepTo(local, tpause);

      // do after-step operations e.g. wait for spike delivery
      synchronizer_->NeuronSyncEnd(local, spikes_lco);
      // if (fmod(t, dt_io) == 0) {  /*output*/ }

      if (has_scheduler) {
        // increment scheduler counter to allow it to look for next job
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
        fprintf(stderr, "~~ %d.%d.%d ~~ before hpx_lco_sema(locality::neurons_scheduler_sema_) 7\n", 
              wrappers::MyRankId(), wrappers::MyThreadId(), wrappers::MyLightWeightThreadId());
#endif
        hpx_lco_sema_v_sync(locality::neurons_scheduler_sema_);
        //hpx_lco_sema_v(locality::neurons_scheduler_sema_, HPX_NULL);
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
        fprintf(stderr, "~~ %d.%d.%d ~~ after hpx_lco_sema(locality::neurons_scheduler_sema_) 7\n", 
              wrappers::MyRankId(), wrappers::MyThreadId(), wrappers::MyLightWeightThreadId());
#endif

        // re-add this job to queue, to be picked up again later
/*
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
        fprintf(stderr, "~~ %d.%d ~~ before hpx_lco_sema(locality::neurons_progress_mutex_) 8\n", 
                wrappers::MyRankId(), wrappers::MyThreadId());
#endif
*/
        hpx_lco_sema_p(locality::neurons_progress_mutex_);
/*
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
        fprintf(stderr, "~~ %d.%d ~~ after hpx_lco_sema(locality::neurons_progress_mutex_) 8\n", 
                wrappers::MyRankId(), wrappers::MyThreadId());
#endif
*/
#ifdef PRINT_TIME_DEPENDENCY
        fprintf(stderr,
                "-- %d.%d -- neuron %d goes back to set with time  %.4f\n",
                wrappers::MyRankId(), wrappers::MyThreadId(),
                local->soma_->gid_, nt->_t);
        TimeDependencySynchronizer::PrintDependencies(local);
#endif

        locality::neurons_progress_->insert(
            std::make_pair(nt->_t, step_trigger));

/*
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
      d fprintf(stderr, "~~ %d.%d ~~ before hpx_lco_sema(locality::neurons_progress_mutex_) 9\n", 
                wrappers::MyRankId(), wrappers::MyThreadId());
#endif
*/
        hpx_lco_sema_v_sync(locality::neurons_progress_mutex_);
/*
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
        fprintf(stderr, "~~ %d.%d ~~ after  hpx_lco_sema(locality::neurons_progress_mutex_) 9\n", 
                wrappers::MyRankId(), wrappers::MyThreadId());
#endif
*/
      }
    }
  }
  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
#if defined(PRINT_TIME_DEPENDENCY_STEP_SIZE) or defined(PRINT_TIME_DEPENDENCY)
  if (has_scheduler)
    // TODO it needs to inform others, to be sure they can step?
    fprintf(stderr, "-- %d.%d.%d -- Neuron %d finished.\n", 
            wrappers::MyRankId(), wrappers::MyThreadId(), wrappers::MyLightWeightThreadId(), local->soma_->gid_);
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
    std::queue<hpx_t> empty_queue;
    std::swap(*neurox::locality::neurons_progress_queue_, empty_queue);
    delete neurox::locality::neurons_progress_queue_;
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
