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

void TimeDependencySynchronizer::BeforeSteps(Branch* b) {
  if (b->soma_) {
    TimeDependencies* time_dependencies =
        (TimeDependencies*)b->soma_->synchronizer_neuron_info_;

    // inform time dependants that must be notified in this step
    time_dependencies->SendSteppingNotification(
        b->nt_->_t, b->nt_->_dt, b->soma_->gid_, b->soma_->synapses_);
    // wait until Im sure I can start and finalize this step at t+dt
    time_dependencies->WaitForTimeDependencyNeurons(b->nt_->_t, b->nt_->_dt,
                                                    b->soma_->gid_);
  }
}

double TimeDependencySynchronizer::GetMaxStep(Branch* branch) {
  // at every step we check for notification intervals
  return branch->nt_->_dt;
}

void TimeDependencySynchronizer::AfterReceiveSpikes(Branch* b, hpx_t target,
                                                    neuron_id_t pre_neuron_id,
                                                    spike_time_t spike_time,
                                                    spike_time_t max_time) {
  // inform soma of this neuron of new time dependency update
  hpx_t top_branch_addr = b->soma_ ? target : b->branch_tree_->top_branch_addr_;
  if (b->soma_) {
    TimeDependencies* time_dependencies =
        (TimeDependencies*)b->soma_->synchronizer_neuron_info_;
    time_dependencies->UpdateTimeDependency(pre_neuron_id, max_time);
  } else {
    TimeDependencySynchronizer::TimeDependencies::SendTimeUpdateMessage(
        top_branch_addr, HPX_NULL, pre_neuron_id, max_time);
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
    spike_time_t maxTimeAllowed =
        t + teps + s->min_delay_ + neuron->refractory_period_;

    /* reset LCO to be used next. any spike or step notification
     * happening after must wait for this spike delivery */
    hpx_lco_wait_reset(s->previous_spike_lco_);

    hpx_action_t add_spike_action = input_params_->locality_comm_reduce_
                                        ? Branch::AddSpikeEventLocality
                                        : Branch::AddSpikeEvent;
    hpx_call(s->synapse_addr_, add_spike_action, s->previous_spike_lco_,
             &neuron->gid_, sizeof(neuron_id_t), &tt, sizeof(spike_time_t),
             &maxTimeAllowed, sizeof(spike_time_t));

#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
    printf(
        "Neuron::sendSpikes: gid %d at time %.3f, informs gid %d of next notif "
        "time =%.3f\n",
        this->gid_, tt, s->destination_gid_, t, s->next_notification_time_);
#endif
  }
  return HPX_NULL;
}

TimeDependencySynchronizer::TimeDependencies::TimeDependencies() {
  libhpx_cond_init(&this->dependencies_wait_condition_);
  libhpx_mutex_init(&this->dependencies_lock_);
  this->dependencies_time_neuron_waits_for_ = 0;  // 0 means not waiting
}

TimeDependencySynchronizer::TimeDependencies::~TimeDependencies() {
  libhpx_cond_destroy(&this->dependencies_wait_condition_);
  libhpx_mutex_destroy(&this->dependencies_lock_);
}

size_t TimeDependencySynchronizer::TimeDependencies::GetDependenciesCount() {
  size_t size = -1;
  libhpx_mutex_lock(&this->dependencies_lock_);
  size = dependencies_map_.size();
  libhpx_mutex_unlock(&this->dependencies_lock_);
  return size;
}

void TimeDependencySynchronizer::TimeDependencies::IncreseDependenciesTime(
    floble_t t) {
  libhpx_mutex_lock(&this->dependencies_lock_);
  for (auto& dependency : dependencies_map_) dependency.second += t;
  libhpx_mutex_unlock(&this->dependencies_lock_);
}

floble_t
TimeDependencySynchronizer::TimeDependencies::GetDependenciesMinTime() {
  // if no dependencies, walk to the end of the simulation
  if (dependencies_map_.empty()) return input_params_->tstop_;

  return std::min_element(dependencies_map_.begin(), dependencies_map_.end(),
                          [](pair<neuron_id_t, floble_t> const& lhs,
                             pair<neuron_id_t, floble_t> const& rhs) {
                            return lhs.second < rhs.second;
                          })
      ->second;
}

void TimeDependencySynchronizer::TimeDependencies::UpdateTimeDependency(
    neuron_id_t src_gid, floble_t dependency_notification_time,
    neuron_id_t my_gid, bool initialization_phase) {
  libhpx_mutex_lock(&this->dependencies_lock_);
  if (initialization_phase) {
    if (dependencies_map_.find(src_gid) ==
        dependencies_map_
            .end())  // in execution without branching this is always true
      dependencies_map_[src_gid] = dependency_notification_time;
    else
      dependencies_map_.at(src_gid) =
          std::max(dependencies_map_.at(src_gid), dependency_notification_time);
  } else {
    assert(dependencies_map_.find(src_gid) != dependencies_map_.end());
    if (dependencies_map_.at(src_gid) <
        dependency_notification_time)  // order of msgs is not guaranteed so
                                       // take
                                       // only last update (highest time value)
    {
      dependencies_map_.at(src_gid) = dependency_notification_time;
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
      printf(
          "-- %d (msg from %d) updates dependenciesMap(%d)=%.11f, notif "
          "time=%.11f, getDependenciesMinTime()=%.11f\n",
          my_gid, src_gid, src_gid, dependencies_map_.at(srcGid),
          dependency_notification_time_, GetDependenciesMinTime());
#endif
    } else {
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
      printf(
          "-- %d (msg from %d) DOES NOT UPDATE dependenciesMap(%d)=%.11f, "
          "notif time=%.11f, getDependenciesMinTime()=%.11f\n",
          my_gid, src_gid, src_gid, dependenciesMap.at(src_gid),
          dependency_notification_time, GetDependenciesMinTime());
#endif
    }

    if (dependencies_time_neuron_waits_for_ >
        0)  // if neuron is waiting for a dependencies time update
      if (GetDependenciesMinTime() >=
          dependencies_time_neuron_waits_for_)  // and new min time allows
                                                // neuron to
                                                // proceed
      {
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
        printf(
            "-- %d (msg from %d) wakes up producer, "
            "GetDependenciesMinTime()=%.11f >= t+dt=%.11f\n",
            my_gid, src_gid, GetDependenciesMinTime(),
            dependencies_time_neuron_waits_for_);
#endif
        dependencies_time_neuron_waits_for_ = 0;  // mark neuron as not asleep
                                                  // anymore
        libhpx_cond_broadcast(
            &this->dependencies_wait_condition_);  // wake up neuron
      }
  }
  libhpx_mutex_unlock(&this->dependencies_lock_);
}

void TimeDependencySynchronizer::TimeDependencies::SendTimeUpdateMessage(
    hpx_t top_branch_addr, hpx_t lco, neuron_id_t preneuron_id,
    spike_time_t max_time, bool init_phase) {
  const hpx_action_t update_time_dep_action =
      input_params_->locality_comm_reduce_
          ? TimeDependencySynchronizer::UpdateTimeDependencyLocality
          : TimeDependencySynchronizer::UpdateTimeDependency;

  hpx_call(top_branch_addr, update_time_dep_action, lco, &preneuron_id,
           sizeof(neuron_id_t), &max_time, sizeof(spike_time_t), &init_phase,
           sizeof(bool));
}
void TimeDependencySynchronizer::TimeDependencies::WaitForTimeDependencyNeurons(
    floble_t t, floble_t dt, int gid) {
  // if I have no dependencies... I'm free to go!
  if (dependencies_map_.size() == 0) return;

#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
  printf("== %d enters TimeDependencies::waitForTimeDependencyNeurons\n", gid);
#endif
  libhpx_mutex_lock(&this->dependencies_lock_);
  if (GetDependenciesMinTime() < t + dt)  // if I cant proceed
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
  assert(GetDependenciesMinTime() >= t + dt);
  libhpx_mutex_unlock(&this->dependencies_lock_);
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
  printf("== %d leaves TimeDependencies::waitForTimeDependencyNeurons\n", gid);
#endif
}

void TimeDependencySynchronizer::TimeDependencies::SendSteppingNotification(
    floble_t t, floble_t dt, int gid, std::vector<Neuron::Synapse*>& synapses) {
  for (Neuron::Synapse*& s : synapses)
    /* if in this time step (-teps to give or take few nanosecs for
     * correction of floating point time roundings) */
    if (s->next_notification_time_ - kTEps <= t + dt) {
      // must have been covered by previous steps
      assert(s->next_notification_time_ >= t);
      s->next_notification_time_ =
          t + s->min_delay_ * TimeDependencies::kNotificationIntervalRatio;
      spike_time_t max_time_allowed =
          t + TimeDependencies::kTEps + s->min_delay_;

      /* wait for previous synapse to be delivered (if any) before telling
       * post-syn neuron to proceed in time */
      hpx_lco_wait(s->previous_spike_lco_);

      TimeDependencySynchronizer::TimeDependencies::SendTimeUpdateMessage(
          s->synapse_soma_addr_, HPX_NULL, gid, max_time_allowed);

#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
      printf("## %d notifies %d he can proceed up to %.6fms\n", gid,
             s->destination_gid_, max_time_allowed);
#endif
    }
}

hpx_action_t TimeDependencySynchronizer::UpdateTimeDependency = 0;
int TimeDependencySynchronizer::UpdateTimeDependency_handler(const int nargs,
                                                             const void* args[],
                                                             const size_t[]) {
  NEUROX_MEM_PIN(Branch);
  assert(nargs == 2 || nargs == 3);

  // auto source = libhpx_parcel_get_source(p);
  const neuron_id_t pre_neuron_id = *(const neuron_id_t*)args[0];
  const spike_time_t max_time = *(const spike_time_t*)args[1];
  const bool init_phase = nargs == 3 ? *(const bool*)args[2] : false;

  assert(local->soma_);
  assert(local->soma_->synchronizer_neuron_info_);
  TimeDependencies* time_dependencies =
      (TimeDependencies*)local->soma_->synchronizer_neuron_info_;
  time_dependencies->UpdateTimeDependency(
      pre_neuron_id, (floble_t)max_time, local->soma_ ? local->soma_->gid_ : -1,
      init_phase);
  NEUROX_MEM_UNPIN
}

hpx_action_t TimeDependencySynchronizer::UpdateTimeDependencyLocality = 0;
int TimeDependencySynchronizer::UpdateTimeDependencyLocality_handler(
    const int nargs, const void* args[], const size_t sizes[]) {
  NEUROX_MEM_PIN(uint64_t);
  assert(nargs == 2 || nargs == 3);
  const neuron_id_t pre_neuron_id = *(const neuron_id_t*)args[0];
  vector<hpx_t>& branch_soma_addrs =
      neurox::locality::netcons_somas_->at(pre_neuron_id);
  hpx_t spikes_lco = hpx_lco_and_new(branch_soma_addrs.size());
  for (hpx_t& soma_addr : branch_soma_addrs)
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
