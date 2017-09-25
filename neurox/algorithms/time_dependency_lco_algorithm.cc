#include "neurox/algorithms/time_dependency_lco_algorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

constexpr floble_t
    TimeDependencyLCOAlgorithm::TimeDependencies::kNotificationIntervalRatio;
constexpr double TimeDependencyLCOAlgorithm::TimeDependencies::kTEps;

TimeDependencyLCOAlgorithm::TimeDependencyLCOAlgorithm() {}

TimeDependencyLCOAlgorithm::~TimeDependencyLCOAlgorithm() {}

const AlgorithmId TimeDependencyLCOAlgorithm::GetId() {
  return AlgorithmId::kBackwardEulerTimeDependencyLCO;
}

const char* TimeDependencyLCOAlgorithm::GetString() {
  return "BackwardEulerTimeDependencyLCO";
}

void TimeDependencyLCOAlgorithm::Init() {
  Algorithm::FixedStepMethodsInit();
  if (input_params_->allreduce_at_locality_)
    throw std::runtime_error(
        "Cant run BackwardEulerTimeDependencyLCO with allReduceAtLocality\n");

  const int allReducesCount = 0;
  hpx_bcast_rsync(AllreduceAlgorithm::AllReducesInfo::SetReductionsPerCommStep,
                  &allReducesCount, sizeof(int));
}

void TimeDependencyLCOAlgorithm::Clear() {}

double TimeDependencyLCOAlgorithm::Launch() {
  int total_steps = Algorithm::GetTotalStepsCount();
  hpx_time_t now = hpx_time_now();
  neurox::wrappers::CallAllNeurons(Branch::BackwardEuler, &total_steps,
                                   sizeof(int));
  double elapsed_time = hpx_time_elapsed_ms(now) / 1e3;
  input::Debugger::RunCoreneuronAndCompareAllBranches();
  return elapsed_time;
}

void TimeDependencyLCOAlgorithm::Run(Branch* b, const void* args) {
  int steps = *(int*)args;

  if (b->soma_) {
    TimeDependencies* time_dependencies =
        (TimeDependencies*)b->soma_->algorithm_metadata_;

    // fixes crash for Algorithm::All when TimeDependency algorithm starts at
    // t=inputParams->tend*2
    // increase notification and dependencies time
    for (Neuron::Synapse*& s : b->soma_->synapses_)
      s->next_notification_time_ += b->nt_->_t;
    time_dependencies->IncreseDependenciesTime(b->nt_->_t);
  }

  for (int step = 0; step < steps; step++) b->BackwardEulerStep();
// Input::Coreneuron::Debugger::stepAfterStepBackwardEuler(local,
// &nrn_threads[this->nt->id], secondorder); //SMP ONLY

#ifndef NDEBUG
  if (b->soma_) printf("-- neuron %d finished\n", b->soma_->gid_);
#endif
}

void TimeDependencyLCOAlgorithm::StepBegin(Branch* b) {
  if (b->soma_) {
    TimeDependencies* time_dependencies =
        (TimeDependencies*)b->soma_->algorithm_metadata_;
    // inform time dependants that must be notified in this step
    time_dependencies->SendSteppingNotification(
        b->nt_->_t, b->nt_->_dt, b->soma_->gid_, b->soma_->synapses_);
    // wait until Im sure I can start and finalize this step at t+dt
    time_dependencies->WaitForTimeDependencyNeurons(b->nt_->_t, b->nt_->_dt,
                                                    b->soma_->gid_);
  }
}

void TimeDependencyLCOAlgorithm::StepEnd(Branch* b, hpx_t) {
  input::Debugger::SingleNeuronStepAndCompare(&nrn_threads[b->nt_->id], b,
                                              input_params_->second_order_);
}

void TimeDependencyLCOAlgorithm::AfterReceiveSpikes(Branch* b, hpx_t target,
                                                    neuron_id_t pre_neuron_id,
                                                    spike_time_t spike_time,
                                                    spike_time_t max_time) {
  // inform soma of this neuron of new time dependency update
  hpx_t top_branch_addr = b->soma_ ? target : b->branch_tree_->top_branch_addr_;
  if (b->soma_) {
    TimeDependencies* time_dependencies =
        (TimeDependencies*)b->soma_->algorithm_metadata_;
    time_dependencies->UpdateTimeDependency(pre_neuron_id, max_time);
  } else
    hpx_call(top_branch_addr, Branch::UpdateTimeDependency, HPX_NULL,
             &pre_neuron_id, sizeof(neuron_id_t), &max_time,
             sizeof(spike_time_t));
}

hpx_t TimeDependencyLCOAlgorithm::SendSpikes(Neuron* neuron, double tt,
                                             double t) {
  const floble_t notification_ratio =
      TimeDependencyLCOAlgorithm::TimeDependencies::kNotificationIntervalRatio;
  const double teps = TimeDependencyLCOAlgorithm::TimeDependencies::kTEps;

  for (Neuron::Synapse*& s : neuron->synapses_) {
    s->next_notification_time_ =
        t + (s->min_delay_ + neuron->refractory_period_) * notification_ratio;
    spike_time_t maxTimeAllowed =
        t + teps + s->min_delay_ + neuron->refractory_period_;

    hpx_lco_wait_reset(s->previous_spike_lco_);  // reset LCO to be used next
    // any spike or step notification happening after must wait for this spike
    // delivery

    hpx_call(s->branch_addr_, Branch::AddSpikeEvent, s->previous_spike_lco_,
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

TimeDependencyLCOAlgorithm::TimeDependencies::TimeDependencies() {
  libhpx_cond_init(&this->dependencies_wait_condition_);
  libhpx_mutex_init(&this->dependencies_lock_);
  this->dependencies_time_neuron_waits_for_ = 0;  // 0 means not waiting
}

TimeDependencyLCOAlgorithm::TimeDependencies::~TimeDependencies() {
  libhpx_cond_destroy(&this->dependencies_wait_condition_);
  libhpx_mutex_destroy(&this->dependencies_lock_);
}

size_t TimeDependencyLCOAlgorithm::TimeDependencies::GetDependenciesCount() {
  size_t size = -1;
  libhpx_mutex_lock(&this->dependencies_lock_);
  size = dependencies_map_.size();
  libhpx_mutex_unlock(&this->dependencies_lock_);
  return size;
}

void TimeDependencyLCOAlgorithm::TimeDependencies::IncreseDependenciesTime(
    floble_t t) {
  libhpx_mutex_lock(&this->dependencies_lock_);
  for (auto& dependency : dependencies_map_) dependency.second += t;
  libhpx_mutex_unlock(&this->dependencies_lock_);
}

floble_t
TimeDependencyLCOAlgorithm::TimeDependencies::GetDependenciesMinTime() {
  assert(dependencies_map_.size() > 0);
  return std::min_element(dependencies_map_.begin(), dependencies_map_.end(),
                          [](pair<neuron_id_t, floble_t> const& lhs,
                             pair<neuron_id_t, floble_t> const& rhs) {
                            return lhs.second < rhs.second;
                          })
      ->second;
}

void TimeDependencyLCOAlgorithm::TimeDependencies::UpdateTimeDependency(
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

void TimeDependencyLCOAlgorithm::TimeDependencies::WaitForTimeDependencyNeurons(
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

void TimeDependencyLCOAlgorithm::TimeDependencies::SendSteppingNotification(
    floble_t t, floble_t dt, int gid, std::vector<Neuron::Synapse*>& synapses) {
  for (Neuron::Synapse*& s : synapses)
    if (s->next_notification_time_ - kTEps <= t + dt)  // if in this time step
    //(-teps to give or take few nanosecs for correction of floating point time
    // roundings)
    {
      assert(s->next_notification_time_ >=
             t);  // must have been covered by previous steps
      s->next_notification_time_ =
          t + s->min_delay_ * TimeDependencies::kNotificationIntervalRatio;
      spike_time_t max_time_allowed =
          t + TimeDependencies::kTEps + s->min_delay_;

      // wait for previous synapse to be delivered (if any) before telling
      // post-syn neuron to proceed in time
      hpx_lco_wait(s->previous_spike_lco_);
      hpx_call(s->top_branch_addr_, Branch::UpdateTimeDependency, HPX_NULL,
               &gid, sizeof(neuron_id_t), &max_time_allowed,
               sizeof(spike_time_t));

#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
      printf("## %d notifies %d he can proceed up to %.6fms\n", gid,
             s->destination_gid_, max_time_allowed);
#endif
    }
}
