#include "neurox/synchronizers/time_dependency_synchronizer.h"

using namespace neurox;
using namespace neurox::synchronizers;

constexpr floble_t
    TimeDependencySynchronizer::TimeDependencies::kNotificationIntervalRatio;
constexpr double TimeDependencySynchronizer::TimeDependencies::kTEps;

TimeDependencySynchronizer::TimeDependencySynchronizer() {
  assert(
      TimeDependencySynchronizer::TimeDependencies::kNotificationIntervalRatio >
          0 &&
      TimeDependencySynchronizer::TimeDependencies::
              kNotificationIntervalRatio <= 1);
}

TimeDependencySynchronizer::~TimeDependencySynchronizer() {}

const SynchronizerIds TimeDependencySynchronizer::GetId() {
  return SynchronizerIds::kTimeDependency;
}

const char* TimeDependencySynchronizer::GetString() {
  return "TimeDependencySynchronizer";
}

void TimeDependencySynchronizer::ClearLocality() {}

void TimeDependencySynchronizer::NeuronSyncEnd(Branch* b, hpx_t) {
  if (!b->soma_) return;
  const bool has_scheduler = b->soma_->scheduler_step_trigger_;
  const bool finished =
      b->nt_->_t > input_params_->tstop_ - TimeDependencies::kTEps;

  if (!has_scheduler) {
#ifdef PRINT_TIME_DEPENDENCY_NEURON_FINISHED
    if (finished) fprintf(stderr, "-- Neuron %d finished.\n", b->soma_->gid_);
#endif
    return;
  }

  TimeDependencies* time_dependencies =
      (TimeDependencies*)b->soma_->synchronizer_neuron_info_;
  // inform Im at end of step (inside waits for any syn that occured)
  time_dependencies->SendSteppingNotification(b);

  // allow scheduler to pick new job
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
  fprintf(stderr,
          "~~ %d.%d.%d ~~ before "
          "hpx_lco_sema(locality::neurons_scheduler_sema_ %d) 7\n",
          wrappers::MyRankId(), wrappers::MyThreadId(),
          wrappers::MyLightWeightThreadId(), locality::scheduler_sema_counter_);
#endif
  // allow scheduler to pick another neuron
  hpx_lco_sema_v_sync(locality::scheduler_neurons_sema_);
  // hpx_lco_sema_v(locality::scheduler_neurons_sema_, HPX_NULL);

#ifdef PRINT_TIME_DEPENDENCY_MUTEX
  fprintf(stderr,
          "~~ %d.%d.%d ~~ after "
          "hpx_lco_sema(locality::neurons_scheduler_sema_ %d) 7\n",
          wrappers::MyRankId(), wrappers::MyThreadId(),
          wrappers::MyLightWeightThreadId(),
          ++locality::scheduler_sema_counter_);
#endif

  if (!finished)
    // wait until min step is allowed
    time_dependencies->WaitForTimeDependencyNeurons(b, kSchedulerMinStep);

  libhpx_mutex_lock(&locality::scheduler_lock_);
  if (finished) {
    locality::scheduler_remaining_neurons_--;
#ifdef PRINT_TIME_DEPEPENCY_NEURON_FINISHED
    fprintf(stderr, "-- Neuron %d finished. Remaining %d.\n", b->soma_->gid_,
            locality::scheduler_remaining_neurons_);
#endif
  } else {
    // add itself to the set of neurons to be run and inform scheduler
    locality::scheduler_neurons_->insert(
        std::make_pair(t, b->soma_->scheduler_step_trigger_));
#ifdef PRINT_TIME_DEPENDENCY
    fprintf(stderr, "-- %d.%d -- neuron %d goes back to set with time  %.4f\n",
            wrappers::MyRankId(), wrappers::MyThreadId(), b->soma_->gid_,
            b->nt_->_t);
    TimeDependencySynchronizer::PrintDependencies(b);
#endif
  }
  // inform scheduler that i finished or Im back in the queue
  libhpx_cond_broadcast(&locality::scheduler_wait_condition_);
  libhpx_mutex_unlock(&locality::scheduler_lock_);
}

double TimeDependencySynchronizer::PrintDependencies(Branch* b) {
  return TimeDependencies::PrintDependencies(b);
}

double TimeDependencySynchronizer::TimeDependencies::PrintDependencies(
    Branch* b) {
  TimeDependencies* td = (TimeDependencies*)b->soma_->synchronizer_neuron_info_;

  double dependencies_min_time = td->GetDependenciesMinTime();
  fprintf(stderr,
          "-- %d.%d -- neuron %d t=%.4f, GetDependenciesMinTime()=%.4f\n",
          wrappers::MyRankId(), wrappers::MyThreadId(), b->soma_->gid_,
          b->nt_->_t, dependencies_min_time);

  for (size_t d = 0; d < td->GetDependenciesCount(); d++) {
    neuron_id_t gid = td->GetDependenciesKeyAtOffset(d);
    floble_t min_delay = td->GetDependencyMinDelay(gid);
    floble_t max_time = td->GetDependencyMaxTimeAllowed(gid);
    printf(
        "   -- pre-syn neuron id %d min delay %.4f allows stepping to %.4f "
        "ms\n",
        wrappers::MyRankId(), wrappers::MyThreadId(), gid, min_delay, max_time);
  }

  size_t syn_count = b->soma_->GetSynapsesCount();
  for (int i = 0; i < syn_count; i++) {
    Neuron::Synapse* s = b->soma_->GetSynapseAtOffset(i);
    printf("  -- post-syn neuron %d, notif time %.4f\n", s->destination_gid_,
           s->next_notification_time_);
  }

  return dependencies_min_time;
}

void TimeDependencySynchronizer::StepSync(Branch* b, const floble_t dt) {
  assert(b->soma_);
  /*if scheduler, send the last step only, and wait for sent spikes.*/
  TimeDependencies* time_dependencies =
      (TimeDependencies*)b->soma_->synchronizer_neuron_info_;
  time_dependencies->SendSteppingNotification(b, dt);
  const bool has_scheduler = b->soma_->scheduler_step_trigger_;
  if (!has_scheduler) {
    /* wait until Im sure I can start and finalize this step at t+dt */
    time_dependencies->WaitForTimeDependencyNeurons(b, dt);
  }
}

double TimeDependencySynchronizer::GetNeuronMaxStep(Branch* b) {
  // at every step we check for notification intervals
  assert(b->soma_);

  const bool has_scheduler = b->soma_->scheduler_step_trigger_;
  if (!has_scheduler) return input_params_->tstop_;

  /* if a scheduler exists, it steps the last neuron until maximum possible
   * time. If it does not exist, this neuron steps freely until tstop */
  TimeDependencies* time_dependencies =
      (TimeDependencies*)b->soma_->synchronizer_neuron_info_;

  // get or wait for a time to step to value, that is higher than current time
  /*
  #ifdef PRINT_TIME_DEPENDENCY_MUTEX
    fprintf(stderr, "~~ %d.%d ~~ neuron %d before
  libhpx_mutex_lock(&time_dependencies->dependencies_lock_) 10\n",
            wrappers::MyRankId(), wrappers::MyThreadId(), b->soma_->gid_);
  #endif
  */
  libhpx_mutex_lock(&time_dependencies->dependencies_lock_);
  double dep_min_time = time_dependencies->GetDependenciesMinTime();
  libhpx_mutex_unlock(&time_dependencies->dependencies_lock_);
  /*
  #ifdef PRINT_TIME_DEPENDENCY_MUTEX
    fprintf(stderr, "~~ %d.%d ~~ neuron %d after
  libhpx_mutex_unlock(&time_dependencies->dependencies_lock_) 11\n",
            wrappers::MyRankId(), wrappers::MyThreadId(), b->soma_->gid_);
  #endif
  */
  // step can be zero if notifications from dependencies didn't arrive yet!
  double step_size = dep_min_time - b->nt_->_t;
  assert(step_size >= kSchedulerMinStep);
  return step_size;
}

void TimeDependencySynchronizer::AfterReceiveSpikes(
    Branch* b, hpx_t target, neuron_id_t pre_neuron_id, spike_time_t,
    spike_time_t dependency_time) {
  // inform soma of this neuron of new time dependency update
  if (b->soma_) {
    TimeDependencies* time_dependencies =
        (TimeDependencies*)b->soma_->synchronizer_neuron_info_;
    time_dependencies->UpdateTimeDependency(pre_neuron_id, dependency_time);
  } else {
    // this branch will tell his some of the time-dependency update
    // Reminder: locality-reduction was already performed to receive this
    // synapse
    hpx_t top_branch_addr =
        b->soma_ ? target : b->branch_tree_->top_branch_addr_;
    hpx_call(top_branch_addr, TimeDependencySynchronizer::UpdateTimeDependency,
             HPX_NULL, &pre_neuron_id, sizeof(neuron_id_t), &dependency_time,
             sizeof(spike_time_t));
  }
}

hpx_t TimeDependencySynchronizer::SendSpikes(Neuron* neuron, double tt,
                                             double t) {
  const floble_t notification_ratio =
      TimeDependencySynchronizer::TimeDependencies::kNotificationIntervalRatio;
  const double teps = TimeDependencySynchronizer::TimeDependencies::kTEps;

  size_t syn_count = neuron->GetSynapsesCount();

  if (input_params_->output_comm_count_) {
    hpx_lco_sema_p(Statistics::CommCount::mutex);
    Statistics::CommCount::point_to_point_count += syn_count;
    hpx_lco_sema_v_sync(Statistics::CommCount::mutex);
  }

  for (int i = 0; i < syn_count; i++) {
    Neuron::Synapse* s = neuron->GetSynapseAtOffset(i);

    /* reminder, s->min_delay_ is the syn. min-delay to branch or locality*/
    s->next_notification_time_ =
        t + (s->min_delay_ + neuron->refractory_period_) * notification_ratio;
    spike_time_t min_time_before_spiking =
        t + teps + neuron->refractory_period_;

#ifdef PRINT_TIME_DEPENDENCY
    fprintf(stderr,
            "-- %d.%d -- neuron %d notifies %d. spikes at time %.4f. next "
            "notif time %.4f\n",
            wrappers::MyRankId(), wrappers::MyThreadId(), neuron->gid_,
            s->destination_gid_, t, s->next_notification_time_);
#endif

    /* reset LCO to be used next. any spike or step notification
     * happening after must wait for this spike delivery */
    hpx_lco_wait_reset(s->previous_synapse_lco_);

    hpx_action_t add_spike_action = input_params_->locality_comm_reduce_
                                        ? Branch::AddSpikeEventLocality
                                        : Branch::AddSpikeEvent;
    hpx_call(s->branch_addr_, add_spike_action, s->previous_synapse_lco_,
             &neuron->gid_, sizeof(neuron_id_t), &tt, sizeof(spike_time_t),
             &min_time_before_spiking, sizeof(spike_time_t));

#ifdef PRINT_TIME_DEPENDENCY
    fprintf(stderr,
            "-- %d.%d -- neuron gid %d spikes at time %.3f, informs gid %d of "
            "next notif "
            "time =%.3f\n",
            wrappers::MyRankId(), wrappers::MyThreadId(), neuron->gid_, tt,
            s->destination_gid_, t, s->next_notification_time_);
#endif
  }
  return HPX_NULL;
}

double TimeDependencySynchronizer::LocalitySyncInterval() {
  // -1 means advance last neuron first
  return input_params_->neurons_scheduler_ ? -1 : 0;
}

TimeDependencySynchronizer::TimeDependencies::TimeDependencies()
    : dependencies_min_delay_linear_(nullptr),
      dependencies_max_time_allowed_linear_(nullptr) {
  libhpx_cond_init(&this->dependencies_wait_condition_);
  libhpx_mutex_init(&this->dependencies_lock_);
  this->dependencies_time_neuron_waits_for_ = 0;  // 0 means not waiting
  this->last_notification_time_ = -1;  // <0 so first step is unconfirmed
}

TimeDependencySynchronizer::TimeDependencies::~TimeDependencies() {
  libhpx_cond_destroy(&this->dependencies_wait_condition_);
  libhpx_mutex_destroy(&this->dependencies_lock_);
}

floble_t TimeDependencySynchronizer::TimeDependencies::GetDependencyMinDelay(
    neuron_id_t gid) {
  return dependencies_min_delay_linear_
             ? *dependencies_min_delay_linear_->At(gid)
             : dependencies_min_delay_.at(gid);
}

void TimeDependencySynchronizer::TimeDependencies::SetDependencyMaxTimeAllowed(
    neuron_id_t gid, floble_t v) {
  if (dependencies_max_time_allowed_linear_)
    *dependencies_max_time_allowed_linear_->At(gid) = v;
  else
    dependencies_max_time_allowed_.at(gid) = v;
}

floble_t
TimeDependencySynchronizer::TimeDependencies::GetDependencyMaxTimeAllowed(
    neuron_id_t gid) {
  return dependencies_max_time_allowed_linear_
             ? *dependencies_max_time_allowed_linear_->At(gid)
             : dependencies_max_time_allowed_.at(gid);
}

neuron_id_t
TimeDependencySynchronizer::TimeDependencies::GetDependenciesKeyAtOffset(
    size_t d) {
  if (dependencies_max_time_allowed_linear_) {
    assert(dependencies_max_time_allowed_linear_->Count() ==
           dependencies_max_time_allowed_linear_->Count());
    return dependencies_max_time_allowed_linear_->KeyAt(d);
  }

  assert(dependencies_max_time_allowed_.size() ==
         dependencies_min_delay_.size());
  return std::next(dependencies_max_time_allowed_.begin(), d)->first;
}

size_t TimeDependencySynchronizer::TimeDependencies::GetDependenciesCount() {
  size_t size = dependencies_min_delay_linear_
                    ? dependencies_min_delay_linear_->Count()
                    : dependencies_min_delay_.size();
  assert(dependencies_min_delay_.size() ==
         dependencies_max_time_allowed_.size());
  return size;
}

floble_t
TimeDependencySynchronizer::TimeDependencies::GetDependenciesMinTime() {
  // if no dependencies, walk to the end of the simulation
  if (GetDependenciesCount() == 0) return input_params_->tstop_;

  if (dependencies_max_time_allowed_linear_) {
    floble_t* max_times = dependencies_max_time_allowed_linear_->Values();
    return *std::min_element(
        max_times,
        max_times + dependencies_max_time_allowed_linear_->ValuesCount());
  }
  return std::min_element(dependencies_max_time_allowed_.begin(),
                          dependencies_max_time_allowed_.end(),
                          [](pair<neuron_id_t, floble_t> const& lhs,
                             pair<neuron_id_t, floble_t> const& rhs) {
                            return lhs.second < rhs.second;
                          })
      ->second;
}

void TimeDependencySynchronizer::TimeDependencies::UpdateTimeDependency(
    neuron_id_t src_gid, floble_t dependency_time, neuron_id_t my_gid,
    bool init_phase) {
  /*
  #ifdef PRINT_TIME_DEPENDENCY_MUTEX
    fprintf(stderr, "~~ %d.%d ~~ before
  libhpx_mutex_lock(&this->dependencies_lock_) 16\n", wrappers::MyRankId(),
  wrappers::MyThreadId()); #endif
  */
  libhpx_mutex_lock(&this->dependencies_lock_);
  /*
  #ifdef PRINT_TIME_DEPENDENCY_MUTEX
    fprintf(stderr, "~~ %d.%d ~~ after
  libhpx_mutex_lock(&this->dependencies_lock_) 16\n", wrappers::MyRankId(),
  wrappers::MyThreadId()); #endif
  */

  /* Reminder: if init_phase is true: dependency_time is the min delay
   * otherwise: dependency_time is the step time or spike time + refraction
   * (time this neuron can advance to = step + refraction + min delay of all
   * synapses)
   */

  if (init_phase) {
    // only non-linear structs are available during data loading
    assert(dependencies_max_time_allowed_linear_ == nullptr);
    assert(dependencies_min_delay_linear_ == nullptr);
    const floble_t max_time_allowed =
        input_params_->tstart_ + dependency_time + kTEps;

    // in execution without branching this is always true
    if (dependencies_max_time_allowed_.find(src_gid) ==
        dependencies_max_time_allowed_.end()) {
      dependencies_max_time_allowed_[src_gid] = max_time_allowed;
      dependencies_min_delay_[src_gid] = dependency_time;
    } else {
      // if pre-syn connects to several branches or instances, take smallest
      // delay and max-time-allowed;
      dependencies_max_time_allowed_.at(src_gid) = std::min(
          dependencies_max_time_allowed_.at(src_gid), max_time_allowed);
      dependencies_min_delay_.at(src_gid) =
          std::min(dependencies_min_delay_.at(src_gid), dependency_time);
    }
  } else {
    const floble_t max_time_allowed =
        dependency_time + GetDependencyMinDelay(src_gid) + kTEps;

    // order of msgs is not guaranteed so take only last update (highest time)
    if (max_time_allowed > GetDependencyMaxTimeAllowed(src_gid)) {
      SetDependencyMaxTimeAllowed(src_gid, max_time_allowed);
      /*
      #ifdef PRINT_TIME_DEPENDENCY
            fprintf(stderr,
                    "-- %d.%d -- %d (msg from %d) updates "
                    "dependencies_max_time_allowed_(%d)=%.11f, notif "
                    "time=%.11f, getDependenciesMinTime()=%.11f\n",
                    wrappers::MyRankId(), wrappers::MyThreadId(),
                    my_gid, src_gid, src_gid,
                    GetDependenciesMaxTimeAlllowed(src_gid),
      max_time_allowed, GetDependenciesMinTime()); #endif
      */
    }

    // if neuron is waiting for a dependencies time update
    if (dependencies_time_neuron_waits_for_ > 0)
      // and new min time allows neuron to proceed
      if (GetDependenciesMinTime() >= dependencies_time_neuron_waits_for_) {
        /*
        #ifdef PRINT_TIME_DEPENDENCY
                fprintf(stderr,
                        "-- %d.%d -- %d (msg from %d) wakes up producer, "
                        "GetDependenciesMinTime()=%.11f >= t+dt=%.11f\n",
                        wrappers::MyRankId(), wrappers::MyThreadId(),
                        my_gid, src_gid, GetDependenciesMinTime(),
                        dependencies_time_neuron_waits_for_);
        #endif
        */
        // mark neuron as not asleep anymore
        dependencies_time_neuron_waits_for_ = 0;
        // wake up neuron
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
        fprintf(
            stderr,
            "~~ %d.%d ~~ neuron %d before "
            "libhpx_cond_broadcast(&this->dependencies_wait_condition_) 17\n",
            wrappers::MyRankId(), wrappers::MyThreadId(), my_gid);
#endif
        libhpx_cond_broadcast(&this->dependencies_wait_condition_);
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
        fprintf(
            stderr,
            "~~ %d.%d ~~ neuron %d after "
            "libhpx_cond_broadcast(&this->dependencies_wait_condition_) 17\n",
            wrappers::MyRankId(), wrappers::MyThreadId(), my_gid);
#endif
      }
  }
  /*
  #ifdef PRINT_TIME_DEPENDENCY_MUTEX
          fprintf(stderr, "~~ %d.%d ~~ neuron %d before
  libhpx_mutex_unlock(&this->dependencies_lock_) 18\n", wrappers::MyRankId(),
  wrappers::MyThreadId(), my_gid); #endif
  */
  libhpx_mutex_unlock(&this->dependencies_lock_);
  /*
  #ifdef PRINT_TIME_DEPENDENCY_MUTEX
          fprintf(stderr, "~~ %d.%d ~~ neuron %d after
  libhpx_mutex_unlock(&this->dependencies_lock_) 18\n", wrappers::MyRankId(),
  wrappers::MyThreadId(), my_gid); #endif
  */
}

void TimeDependencySynchronizer::TimeDependencies::WaitForTimeDependencyNeurons(
    Branch* b, const floble_t dt) {
  // if neuron has no dependencies... no need to wait
  if (GetDependenciesCount() == 0) return;

  const floble_t t = b->nt_->_t;

  // if this is an "end of execution notification"... no need to wait
  if (t > input_params_->tstop_ - TimeDependencies::kTEps) return;

#ifdef PRINT_TIME_DEPENDENCY_MUTEX
  fprintf(
      stderr,
      "~~ %d.%d ~~ before libhpx_mutex_lock(&this->dependencies_lock_) 19\n",
      wrappers::MyRankId(), wrappers::MyThreadId());
#endif
  libhpx_mutex_lock(&this->dependencies_lock_);
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
  fprintf(stderr,
          "~~ %d.%d ~~ after libhpx_mutex_lock(&this->dependencies_lock_) 19\n",
          wrappers::MyRankId(), wrappers::MyThreadId());
#endif
  if (GetDependenciesMinTime() + kTEps < t + dt)  // if I cant proceed
  {
#ifdef PRINT_TIME_DEPENDENCY
    fprintf(stderr,
            "== %d.%d == %d cant proceed and sleeps: "
            "GetDependenciesMinTime()=%.11f < "
            "t+dt=%.11f\n",
            wrappers::MyRankId(), wrappers::MyThreadId(), b->soma_->gid_,
            GetDependenciesMinTime(), t + dt);
    PrintDependencies(b);
#endif
    // mark this neuron as asleep waiting for a given min dependencies time
    dependencies_time_neuron_waits_for_ = t + dt;
    // release dependenciesLock and sleep until woken up
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
    fprintf(stderr,
            "~~ %d.%d ~~ before "
            "libhpx_cond_wait(&this->dependencies_wait_condition_, "
            "&this->dependencies_lock_); 20\n",
            wrappers::MyRankId(), wrappers::MyThreadId());
#endif

    // wait until my dependencies allow me to step
    libhpx_cond_wait(&this->dependencies_wait_condition_,
                     &this->dependencies_lock_);

#ifdef PRINT_TIME_DEPENDENCY_MUTEX
    fprintf(stderr,
            "~~ %d.%d ~~ after "
            "libhpx_cond_wait(&this->dependencies_wait_condition_, "
            "&this->dependencies_lock_); 20\n",
            wrappers::MyRankId(), wrappers::MyThreadId());
#endif

#ifdef PRINT_TIME_DEPENDENCY
    fprintf(stderr, "== %d.%d == %d wakes up: getDependenciesMinTime()=%.11f\n",
            wrappers::MyRankId(), wrappers::MyThreadId(), b->soma_->gid_,
            GetDependenciesMinTime());
#endif
  } else {
#ifdef PRINT_TIME_DEPENDENCY
    fprintf(stderr,
            "== %d.%d == %d proceeds: getDependenciesMinTime()=%.11f >= "
            "t+dt=%.11f\n",
            wrappers::MyRankId(), wrappers::MyThreadId(), b->soma_->gid_,
            GetDependenciesMinTime(), t + dt);
#endif
  }
  assert(GetDependenciesMinTime() + kTEps >= t + dt);
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
  fprintf(
      stderr,
      "~~ %d.%d ~~ before libhpx_mutex_unlock(&this->dependencies_lock_) 21\n",
      wrappers::MyRankId(), wrappers::MyThreadId());
#endif
  libhpx_mutex_unlock(&this->dependencies_lock_);
#ifdef PRINT_TIME_DEPENDENCY_MUTEX
  fprintf(
      stderr,
      "~~ %d.%d ~~ after libhpx_mutex_unlock(&this->dependencies_lock_) 21\n",
      wrappers::MyRankId(), wrappers::MyThreadId());
#endif
}

void TimeDependencySynchronizer::TimeDependencies::SendSteppingNotification(
    Branch* b) {
  return SendSteppingNotification(b, -1);  //-1 means send anyway
}

void TimeDependencySynchronizer::TimeDependencies::SendSteppingNotification(
    Branch* b, const floble_t dt) {
  if (b->soma_->GetSynapsesCount() == 0) return;

  const floble_t t = b->nt_->_t;

  // avoid sending repeated notifications
  // (useful on var dt or fixed dt with synchronizer when doest step)
  if (t == this->last_notification_time_) return;

  Neuron* neuron = b->soma_;
  const neuron_id_t gid = b->soma_->gid_;

#ifdef PRINT_TIME_DEPENDENCY
  fprintf(stderr, "-- %d.%d.%d -- neuron %d before sending spikes\n",
          wrappers::MyRankId(), wrappers::MyThreadId(),
          wrappers::MyLightWeightThreadId(), gid);
#endif

  // TODO: avoid sending messages that are too close
  // if (t - this->last_notification_time_ < neurox::min_synaptic_delay/2)
  // return;

  this->last_notification_time_ = t;

  const hpx_action_t update_time_dep_action =
      input_params_->locality_comm_reduce_
          ? TimeDependencySynchronizer::UpdateTimeDependencyLocality
          : TimeDependencySynchronizer::UpdateTimeDependency;

  size_t syn_count = neuron->GetSynapsesCount();
  if (input_params_->output_comm_count_) {
    hpx_lco_sema_p(Statistics::CommCount::mutex);
    Statistics::CommCount::point_to_point_count += syn_count;
    hpx_lco_sema_v_sync(Statistics::CommCount::mutex);
  }

  for (int i = 0; i < syn_count; i++) {
    Neuron::Synapse* s = neuron->GetSynapseAtOffset(i);

    /* if in this time step (-teps to give or take few nanosecs for
     * correction of floating point time roundings) ;
     * -1 is a flag, means send anyway */
    if (dt == -1 || s->next_notification_time_ - kTEps <= t + dt) {
      s->next_notification_time_ =
          t + s->min_delay_ * TimeDependencies::kNotificationIntervalRatio;

      // commented: for variable dt, one can jump ahead of notification time
      // assert(s->next_notification_time_ >= t);

      // Wait for previous synapse to be delivered, if any (does not reset)
      hpx_lco_wait(s->previous_synapse_lco_);
      hpx_call(s->soma_or_locality_addr_, update_time_dep_action, HPX_NULL,
               &gid, sizeof(neuron_id_t), &t, sizeof(spike_time_t));
    }
  }
  /*
  #ifdef PRINT_TIME_DEPENDENCY
        fprintf(stderr, "-- %d.%d -- neuron %d after sending
  spikes/notification)\n", wrappers::MyRankId(), wrappers::MyThreadId(),
  gid); #endif
  */
}

hpx_action_t TimeDependencySynchronizer::UpdateTimeDependency = 0;
int TimeDependencySynchronizer::UpdateTimeDependency_handler(const int nargs,
                                                             const void* args[],
                                                             const size_t[]) {
  NEUROX_MEM_PIN(Branch);
  assert(nargs == 2 || nargs == 3);

  // if this neurons already finished, discard dependency step notifications
  if (local->nt_->_t > input_params_->tstop_ - TimeDependencies::kTEps)
    NEUROX_MEM_UNPIN;

  const neuron_id_t pre_neuron_id = *(const neuron_id_t*)args[0];
  const spike_time_t dependency_time = *(const spike_time_t*)args[1];
  const bool init_phase = nargs == 3 ? *(const bool*)args[2] : false;

  /*
  #ifdef PRINT_TIME_DEPENDENCY
    fprintf(stderr, "-- %d.%d -- neuron %d is notified by %d of time %.4f\n",
            wrappers::MyRankId(), wrappers::MyThreadId(),
            local->soma_->gid_, pre_neuron_id, dependency_time);
    assert(local->soma_ && local->soma_->synchronizer_neuron_info_);
  #endif
  */
  TimeDependencies* time_dependencies =
      (TimeDependencies*)local->soma_->synchronizer_neuron_info_;
  time_dependencies->UpdateTimeDependency(
      pre_neuron_id, (floble_t)dependency_time,
      local->soma_ ? local->soma_->gid_ : -1, init_phase);
  NEUROX_MEM_UNPIN
}

hpx_action_t TimeDependencySynchronizer::UpdateTimeDependencyLocality = 0;
int TimeDependencySynchronizer::UpdateTimeDependencyLocality_handler(
    const int nargs, const void* args[], const size_t sizes[]) {
  NEUROX_MEM_PIN(uint64_t);
  assert(nargs == 2 || nargs == 3);
  const neuron_id_t pre_neuron_id = *(const neuron_id_t*)args[0];
  vector<hpx_t>& somas_addrs =
      neurox::locality::netcons_somas_->at(pre_neuron_id);
  hpx_t spikes_lco = hpx_lco_and_new(somas_addrs.size());

  // inform all somas that pre_neuron connects to about the pre-neuron step
  for (hpx_t& soma_addr : somas_addrs)
    if (nargs == 2)
      hpx_call(soma_addr, TimeDependencySynchronizer::UpdateTimeDependency,
               spikes_lco, args[0], sizes[0], args[1], sizes[1]);
    else
      hpx_call(soma_addr, TimeDependencySynchronizer::UpdateTimeDependency,
               spikes_lco, args[0], sizes[0], args[1], sizes[1], args[2],
               sizes[2]);
  hpx_lco_wait(spikes_lco);
  hpx_lco_delete(spikes_lco, HPX_NULL);
  NEUROX_MEM_UNPIN
}

void TimeDependencySynchronizer::RegisterHpxActions() {
  wrappers::RegisterMultipleVarAction(UpdateTimeDependency,
                                      UpdateTimeDependency_handler);
  wrappers::RegisterMultipleVarAction(UpdateTimeDependencyLocality,
                                      UpdateTimeDependencyLocality_handler);
}
