#pragma once
#include "neurox.h"

#include "libhpx/libhpx.h"

using namespace neurox;

namespace neurox {

namespace synchronizers {

class TimeDependencySynchronizer : public Synchronizer {
 public:
  TimeDependencySynchronizer();
  ~TimeDependencySynchronizer();

  const SynchronizerIds GetId() override;
  const char* GetString() override;

  void InitNeuron(Branch*) override;
  void ClearLocality() override;

  double GetNeuronMaxStep(Branch*) override;
  hpx_t SendSpikes(Neuron* b, double tt, double t) override;
  double LocalitySyncInterval();

  void StepBegin(Branch*) override;

  void AfterReceiveSpikes(Branch* local, hpx_t target,
                          neuron_id_t pre_neuron_id, spike_time_t spike_time,
                          spike_time_t max_time) override;

  static hpx_action_t UpdateTimeDependency;
  static hpx_action_t UpdateTimeDependencyLocality;

  static void RegisterHpxActions();  ///> Register all HPX actions

  /// controls time-dependencies from incoming neuron connections
  class TimeDependencies : public SynchronizerNeuronInfo {
   public:
    TimeDependencies();
    ~TimeDependencies();

    void WaitForTimeDependencyNeurons(Branch* b);

    /// inform my outgoing-connection neurons that I stepped
    void SendSteppingNotification(Branch* b);

    /// update time of a given dependency
    void UpdateTimeDependency(neuron_id_t src_gid,
                              floble_t dependency_notification_time,
                              neuron_id_t my_gid = -1,
                              bool initialization = false);

    static void SendTimeUpdateMessage(hpx_t soma_or_locality_addr, hpx_t lco,
                                      neuron_id_t preneuron_id,
                                      spike_time_t max_time,
                                      bool init_phase = false);

    /// get smallest time across all dependencies
    floble_t GetDependenciesMinTime();

    /// get number of dependencies
    size_t GetDependenciesCount();

    /// increase time of all dependencies by 't'
    void IncreseDependenciesTime(floble_t t);

    /// ratio of notification interval (0,1]
    static constexpr const floble_t kNotificationIntervalRatio = 1;

    /// time-epsilon to correct wrong delivery of events due to floating point
    /// rounding
    static constexpr const double kTEps = 1e-8;

   private:
    ///> map of actual time per dependency if
    std::map<neuron_id_t, floble_t> dependencies_map_;
    libhpx_cond_t dependencies_wait_condition_;
    libhpx_mutex_t dependencies_lock_;
    floble_t dependencies_time_neuron_waits_for_;
  };

 private:
  static int UpdateTimeDependency_handler(const int, const void* [],
                                          const size_t[]);
  static int UpdateTimeDependencyLocality_handler(const int, const void* [],
                                                  const size_t[]);
};

};  // synchronizer

};  // neurox
