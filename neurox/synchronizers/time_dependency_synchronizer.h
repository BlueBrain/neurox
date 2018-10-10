#pragma once
#include "neurox.h"

using namespace neurox;

namespace neurox {

namespace synchronizers {

class TimeDependencySynchronizer : public Synchronizer {
 public:
  TimeDependencySynchronizer();
  ~TimeDependencySynchronizer();

  const SynchronizerIds GetId() override;
  const char* GetString() override;

  void NeuronSyncEnd(Branch*, hpx_t) override;
  void ClearLocality() override;

  double GetNeuronMaxStep(Branch*) override;
  hpx_t SendSpikes(Neuron* b, double tt, double t) override;
  double LocalitySyncInterval();

  void StepSync(Branch*, const floble_t dt) override;

  void AfterReceiveSpikes(Branch* local, hpx_t target,
                          neuron_id_t pre_neuron_id, spike_time_t spike_time,
                          spike_time_t dependency_time) override;

  /// minimum step allows by scheduler
  static constexpr const floble_t kSchedulerMinStep = 0.1;

  //// For debugging purposes: outputs dependencies of given branch
  static double PrintDependencies(Branch*);

  static hpx_action_t UpdateTimeDependency;
  static hpx_action_t UpdateTimeDependencyLocality;

  static void RegisterHpxActions();  ///> Register all HPX actions

  /// controls time-dependencies from incoming neuron connections
  class TimeDependencies : public SynchronizerNeuronInfo {
   public:
    TimeDependencies();
    ~TimeDependencies();

    void WaitForTimeDependencyNeurons(Branch* b, const floble_t dt);

    /// inform my outgoing-connection neurons that I stepped
    void SendSteppingNotification(Branch* b);
    void SendSteppingNotification(Branch* b, const floble_t dt);

    /// update time of a given dependency
    void UpdateTimeDependency(neuron_id_t src_gid, floble_t dependency_time,
                              neuron_id_t my_gid = -1, bool init = false);

    /// get dependency min delay
    inline floble_t GetDependencyMinDelay(neuron_id_t gid);

    /// get max time allowed by dependencies
    inline floble_t GetDependencyMaxTimeAllowed(neuron_id_t gid);
    inline void SetDependencyMaxTimeAllowed(neuron_id_t gid, floble_t v);

    /// get smallest time across all dependencies
    floble_t GetDependenciesMinTime();

    /// get number of dependencies
    size_t GetDependenciesCount();
    inline neuron_id_t GetDependenciesKeyAtOffset(size_t d);

    //// For debugging purposes: outputs dependencies of given branch
    static double PrintDependencies(Branch*);

    /// ratio of notification interval (0,1]
    static constexpr const floble_t kNotificationIntervalRatio = 1;

    /// time-epsilon to correct wrong delivery of events due to floating point
    /// rounding
    static constexpr const double kTEps = 1e-8;

    /// controls sleep and waking of neurons after dependencies time-update
    libhpx_mutex_t dependencies_lock_;

    ///> map of synaptic delay per pre-synaptic id
    std::map<neuron_id_t, floble_t> dependencies_min_delay_;
    linear::Map<neuron_id_t, floble_t>* dependencies_min_delay_linear_;

    ///> map of actual time per dependency if
    std::map<neuron_id_t, floble_t> dependencies_max_time_allowed_;
    linear::Map<neuron_id_t, floble_t>* dependencies_max_time_allowed_linear_;

   private:
    /// wait condition that wakes dependencies_lock_
    libhpx_cond_t dependencies_wait_condition_;

    /// time that this neuron waits for, before waking up and continuing
    floble_t dependencies_time_neuron_waits_for_;

    /// time of the last step that is confirmed (useful for var-dt, to avoid
    /// repetition of notification messages
    floble_t last_notification_time_;
  };

 private:
  static int UpdateTimeDependency_handler(const int, const void* [],
                                          const size_t[]);
  static int UpdateTimeDependencyLocality_handler(const int, const void* [],
                                                  const size_t[]);
};

};  // namespace synchronizers

};  // namespace neurox
