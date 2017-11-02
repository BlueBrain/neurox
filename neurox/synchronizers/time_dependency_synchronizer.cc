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

void TimeDependencySynchronizer::InitNeuron(Branch* b) {
  if (b->soma_) {
    TimeDependencies* time_dependencies =
        (TimeDependencies*)b->soma_->synchronizer_neuron_info_;

    /* fixes crash for Synchronizer::All when TimeDependency
     * synchronizer starts at t=inputParams->tend*2 increase
     * notification and dependencies time */
    if (input_params_->synchronizer_ == SynchronizerIds::kBenchmarkAll) {
      for (Neuron::Synapse*& s : b->soma_->synapses_)
        s->next_notification_time_ += b->nt_->_t;
      time_dependencies->IncreseDependenciesTime(b->nt_->_t);
    }
  }
}

void TimeDependencySynchronizer::StepSync(Branch* b, const floble_t dt) {
  assert(b->soma_);
  TimeDependencies* time_dependencies =
      (TimeDependencies*)b->soma_->synchronizer_neuron_info_;

  // inform time dependants that must be notified in this step
  time_dependencies->SendSteppingNotification(b, dt);

  /* wait until Im sure I can start and finalize this step at t+dt:
   * if a scheduler exists, it only allows me to step as far as dependencies
   * allows, so no need to wait for dependencies -- see GetNeuronMaxStep()*/
  const bool has_scheduler = b->soma_->synchronizer_step_trigger_;
  if (!has_scheduler) time_dependencies->WaitForTimeDependencyNeurons(b);
}

double TimeDependencySynchronizer::GetNeuronMaxStep(Branch* b) {
  // at every step we check for notification intervals
  assert(b->soma_);

  const bool has_scheduler = b->soma_->synchronizer_step_trigger_;
  if (!has_scheduler) return input_params_->tstop_;

  /* if a scheduler exists, it steps the last neuron until maximum possible
   * time. If it does not exist, this neuron steps freely until tstop and
   * independently waits for dependencies at every step, via method
   * WaitForTimeDependencyNeurons*/
  TimeDependencies* time_dependencies =
      (TimeDependencies*)b->soma_->synchronizer_neuron_info_;

  // get or wait for a time to step to value, that is higher than current time
  libhpx_mutex_lock(&time_dependencies->dependencies_lock_);
  double dep_min_time = time_dependencies->GetDependenciesMinTime();
  libhpx_mutex_unlock(&time_dependencies->dependencies_lock_);
  double step_size = dep_min_time - b->nt_->_t;

  const double min_step_size = 0.025;
  while (step_size<min_step_size) //wait for dependencies if no step or step too small
  {
    // step notification msg was not processed yet, so this branch will
    // (notify scheduler to pick other neuron,) sleep and use idle cores
    // to process arrival of stepping notifications
    // hpx_lco_sema_p(locality::neurons_progress_mutex_);
    time_dependencies->WaitForTimeDependencyNeurons(b, min_step_size);
    libhpx_mutex_lock(&time_dependencies->dependencies_lock_);
    dep_min_time = time_dependencies->GetDependenciesMinTime();
    libhpx_mutex_unlock(&time_dependencies->dependencies_lock_);
    step_size = dep_min_time - b->nt_->_t;
    // hpx_lco_sema_v_sync(locality::neurons_progress_mutex_);
  }
  assert(step_size>0);
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

  for (Neuron::Synapse*& s : neuron->synapses_) {
    /* reminder, s->min_delay_ is the syn. min-delay to branch or locality*/
    s->next_notification_time_ =
        t + (s->min_delay_ + neuron->refractory_period_) * notification_ratio;
    spike_time_t min_time_before_spiking =
        t + teps + neuron->refractory_period_;

    /* reset LCO to be used next. any spike or step notification
     * happening after must wait for this spike delivery */
    hpx_lco_wait_reset(s->previous_spike_lco_);

    hpx_action_t add_spike_action = input_params_->locality_comm_reduce_
                                        ? Branch::AddSpikeEventLocality
                                        : Branch::AddSpikeEvent;
    hpx_call(s->branch_addr_, add_spike_action, s->previous_spike_lco_,
             &neuron->gid_, sizeof(neuron_id_t), &tt, sizeof(spike_time_t),
             &min_time_before_spiking, sizeof(spike_time_t));

#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
    printf(
        "Neuron::sendSpikes: gid %d at time %.3f, informs gid %d of next notif "
        "time =%.3f\n",
        this->gid_, tt, s->destination_gid_, t, s->next_notification_time_);
#endif
  }
  return HPX_NULL;
}

double TimeDependencySynchronizer::LocalitySyncInterval() {
  return -1;  // means advance last neuron first (see synchronizer.h)
}

TimeDependencySynchronizer::TimeDependencies::TimeDependencies() {
  libhpx_cond_init(&this->dependencies_wait_condition_);
  libhpx_mutex_init(&this->dependencies_lock_);
  this->dependencies_time_neuron_waits_for_ = 0;  // 0 means not waiting
  this->last_notification_time_ = -1; // <0 so first step is unconfirmed
}

TimeDependencySynchronizer::TimeDependencies::~TimeDependencies() {
  libhpx_cond_destroy(&this->dependencies_wait_condition_);
  libhpx_mutex_destroy(&this->dependencies_lock_);
}

size_t TimeDependencySynchronizer::TimeDependencies::GetDependenciesCount() {
  size_t size = -1;
  libhpx_mutex_lock(&this->dependencies_lock_);
  size = dependencies_max_time_allowed_.size();
  libhpx_mutex_unlock(&this->dependencies_lock_);
  return size;
}

void TimeDependencySynchronizer::TimeDependencies::IncreseDependenciesTime(
    floble_t t) {
  libhpx_mutex_lock(&this->dependencies_lock_);
  for (auto& dependency : dependencies_max_time_allowed_)
    dependency.second += t;
  libhpx_mutex_unlock(&this->dependencies_lock_);
}

floble_t
TimeDependencySynchronizer::TimeDependencies::GetDependenciesMinTime() {
  // if no dependencies, walk to the end of the simulation
  // TODO this only works if This is the first sync in the benchmarks (tstop?)
  if (dependencies_max_time_allowed_.empty()) return input_params_->tstop_;

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
  libhpx_mutex_lock(&this->dependencies_lock_);

  /* Reminder: if init_phase is true: dependency_time is the min delay
   * otherwise: dependency_time is the step time or spike time + refraction
   * (time this neuron can advance to = step + refraction + min delay of all
   * synapses)
   */

  if (init_phase) {
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
        dependency_time + dependencies_min_delay_.at(src_gid) + kTEps;

    assert(dependencies_max_time_allowed_.find(src_gid) !=
           dependencies_max_time_allowed_.end());
    // order of msgs is not guaranteed so take only last update (highest time)
    if (dependencies_max_time_allowed_.at(src_gid) < max_time_allowed) {
      dependencies_max_time_allowed_.at(src_gid) = max_time_allowed;
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
      printf(
          "-- %d (msg from %d) updates dependenciesMap(%d)=%.11f, notif "
          "time=%.11f, getDependenciesMinTime()=%.11f\n",
          my_gid, src_gid, src_gid, dependencies_map_.at(srcGid),
          max_time_allowed, GetDependenciesMinTime());
#endif
    } else {
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
      printf(
          "-- %d (msg from %d) DOES NOT UPDATE dependenciesMap(%d)=%.11f, "
          "notif time=%.11f, getDependenciesMinTime()=%.11f\n",
          my_gid, src_gid, src_gid, dependenciesMap.at(src_gid),
          max_time_allowed, GetDependenciesMinTime());
#endif
    }

    // if neuron is waiting for a dependencies time update
    if (dependencies_time_neuron_waits_for_ > 0)
      // and new min time allows neuron to proceed
      if (GetDependenciesMinTime() >= dependencies_time_neuron_waits_for_) {
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
        printf(
            "-- %d (msg from %d) wakes up producer, "
            "GetDependenciesMinTime()=%.11f >= t+dt=%.11f\n",
            my_gid, src_gid, GetDependenciesMinTime(),
            dependencies_time_neuron_waits_for_);
#endif
        // mark neuron as not asleep anymore
        dependencies_time_neuron_waits_for_ = 0;
        // wake up neuron
        libhpx_cond_broadcast(&this->dependencies_wait_condition_);
      }
  }
  libhpx_mutex_unlock(&this->dependencies_lock_);
}

void TimeDependencySynchronizer::TimeDependencies::WaitForTimeDependencyNeurons(
    Branch *b)
{
    return SendSteppingNotification(b, b->nt_->_dt);
}

void TimeDependencySynchronizer::TimeDependencies::WaitForTimeDependencyNeurons(
    Branch *b, const floble_t dt)
{
  // if neuron has no dependencies... no need to wait
  if (dependencies_max_time_allowed_.size() == 0) return;

  // if this is an "end of execution notification"... no need to wait
  const floble_t t = b->nt_->_t;
  if (fabs (t - input_params_->tstop_) < 0.0001 ) return;

#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
  const int gid = b->soma_->gid_;
  printf("== %d enters TimeDependencies::waitForTimeDependencyNeurons\n", gid);
#endif
  libhpx_mutex_lock(&this->dependencies_lock_);
  if (GetDependenciesMinTime() + kTEps < t + dt)  // if I cant proceed
  {
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
    printf(
        "== %d cant proceed and sleeps: GetDependenciesMinTime()=%.11f < "
        "t+dt=%.11f\n",
        gid, GetDependenciesMinTime(), t + dt);
#endif
    // mark this neuron as asleep waiting for a given min dependencies time
    dependencies_time_neuron_waits_for_ = t + dt;
    // release dependenciesLock and sleep until woken up by
    // TimeDependencies::dependenciesWaitCondition
    libhpx_cond_wait(&this->dependencies_wait_condition_,
                     &this->dependencies_lock_);
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
    printf("== %d wakes up: getDependenciesMinTime()=%.11f\n", gid,
           GetDependenciesMinTime());
#endif
  }
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
  else
    printf("== %d proceeds: getDependenciesMinTime()=%.11f >= t+dt=%.11f\n",
           gid, GetDependenciesMinTime(), t + dt);
#endif
  assert(GetDependenciesMinTime() + kTEps >= t + dt);
  libhpx_mutex_unlock(&this->dependencies_lock_);
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
  printf("== %d leaves TimeDependencies::waitForTimeDependencyNeurons\n", gid);
#endif
}

void TimeDependencySynchronizer::TimeDependencies::SendSteppingNotification(
    Branch *b)
{
    return SendSteppingNotification(b, b->nt_->_dt);
}

void TimeDependencySynchronizer::TimeDependencies::SendSteppingNotification(
    Branch* b, const floble_t dt) {
  std::vector<Neuron::Synapse*>& synapses = b->soma_->synapses_;
  if (synapses.empty()) return;

  const floble_t t = b->nt_->_t;
  const neuron_id_t gid = b->soma_->gid_;

  //avoid sending repeated notifications (useful on variable dt)
  if (t == this->last_notification_time_) return;

  //TODO: avoid sending messages that are too close
  //if (t - this->last_notification_time_ < neurox::min_synaptic_delay/2) return;

  this->last_notification_time_ = t;

  const hpx_action_t update_time_dep_action =
      input_params_->locality_comm_reduce_
          ? TimeDependencySynchronizer::UpdateTimeDependencyLocality
          : TimeDependencySynchronizer::UpdateTimeDependency;

  //printf("=============== Neuron %d informs step at %.5f :", gid, t);

  for (Neuron::Synapse*& s : synapses) {
    /* if in this time step (-teps to give or take few nanosecs for
     * correction of floating point time roundings) */
    //printf("%d (", s->destination_gid_);
    if (dt == 0 || s->next_notification_time_ - kTEps <= t + dt) {

      //printf("%f),", s->next_notification_time_);
      s->next_notification_time_ =
        t + s->min_delay_ * TimeDependencies::kNotificationIntervalRatio;


      //commented: for variable dt, one can jump ahead of notification time
      //assert(s->next_notification_time_ >= t);

      /* wait for previous synapse to be delivered (if any) before telling
       * post-syn neuron that I've reached time 't' */
      hpx_lco_wait(s->previous_spike_lco_);

      hpx_call(s->soma_or_locality_addr_, update_time_dep_action, HPX_NULL,
               &gid, sizeof(neuron_id_t), &t, sizeof(spike_time_t));

#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
      printf("## %d notifies %d he that he steppted at %.6fms\n", gid,
             s->destination_gid_, t);
#endif
    }
  }

  //printf("\n");
}

hpx_action_t TimeDependencySynchronizer::UpdateTimeDependency = 0;
int TimeDependencySynchronizer::UpdateTimeDependency_handler(const int nargs,
                                                             const void* args[],
                                                             const size_t[]) {
  NEUROX_MEM_PIN(Branch);
  assert(nargs == 2 || nargs == 3);

  // auto source = libhpx_parcel_get_source(p);
  const neuron_id_t pre_neuron_id = *(const neuron_id_t*)args[0];
  const spike_time_t dependency_time = *(const spike_time_t*)args[1];
  const bool init_phase = nargs == 3 ? *(const bool*)args[2] : false;

  assert(local->soma_ && local->soma_->synchronizer_neuron_info_);
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
