#pragma once
#include "neurox.h"

#include "libhpx/libhpx.h"

using namespace neurox;

namespace neurox {

namespace algorithms {

class TimeDependencyLCOAlgorithm : public Algorithm {
 public:
  TimeDependencyLCOAlgorithm();
  ~TimeDependencyLCOAlgorithm();

  const SyncAlgorithms GetId() override;
  const char* GetString() override;

  void Init() override;
  void Clear() override;
  double Launch() override;

  void StepBegin(Branch*) override;
  void StepEnd(Branch*, hpx_t) override;
  void Run(Branch*, const void*) override;
  hpx_t SendSpikes(Neuron* b, double tt, double t) override;
  void AfterReceiveSpikes(Branch* local, hpx_t target,
                          neuron_id_t pre_neuron_id, spike_time_t spike_time,
                          spike_time_t max_time) override;

  /// controls time-dependencies from incoming neuron connections
  class TimeDependencies : public AlgorithmMetadata {
   public:
    TimeDependencies();
    ~TimeDependencies();

    void WaitForTimeDependencyNeurons(floble_t t, floble_t dt, int gid);

    /// inform my outgoing-connection neurons that I stepped
    void SendSteppingNotification(floble_t t, floble_t dt, int gid,
                                  std::vector<Neuron::Synapse*>& synapses);

    /// update time of a given dependency
    void UpdateTimeDependency(neuron_id_t src_gid,
                              floble_t dependency_notification_time,
                              neuron_id_t my_gid = -1,
                              bool initialization = false);

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
};

};  // algorithm

};  // neurox
