#pragma once
#include "neurox.h"

using namespace neurox;

namespace neurox {

namespace algorithms {

class TimeDependencyLCOAlgorithm : public Algorithm {
 public:
  TimeDependencyLCOAlgorithm();
  ~TimeDependencyLCOAlgorithm();

  const AlgorithmType getType() override;
  const char* getTypeString() override;

  void Init() override;
  void Clear() override;
  double Launch() override;

  void StepBegin(Branch*) override;
  void StepEnd(Branch*, hpx_t) override;
  void Run(Branch*, const void*) override;
  hpx_t SendSpikes(Neuron* b, double tt, double t) override;
  void AfterReceiveSpikes(Branch* local, hpx_t target, neuron_id_t preNeuronId,
                          spike_time_t spikeTime,
                          spike_time_t maxTime) override;

  /// controls time-dependencies from incoming neuron connections
  class TimeDependencies : public AlgorithmMetaData {
   public:
    TimeDependencies();
    ~TimeDependencies();

    void WaitForTimeDependencyNeurons(floble_t t, floble_t dt, int gid);

    /// inform my outgoing-connection neurons that I stepped
    void SendSteppingNotification(floble_t t, floble_t dt, int gid,
                                  std::vector<Neuron::Synapse*>& synapses);

    /// update time of a given dependency
    void UpdateTimeDependency(neuron_id_t srcGid,
                              floble_t dependencyNotificationTime,
                              neuron_id_t myGid = -1,
                              bool initialization = false);

    /// get smallest time across all dependencies
    floble_t GetDependenciesMinTime();

    /// get number of dependencies
    size_t GetDependenciesCount();

    /// increase time of all dependencies by 't'
    void IncreseDependenciesTime(floble_t t);

    /// ratio of notification interval (0,1]
    static floble_t notificationIntervalRatio;

    /// time-epsilon to correct wrong delivery of events due to floating point
    /// rounding
    static double teps;

   private:
    ///> map of actual time per dependency if
    std::map<neuron_id_t, floble_t> dependenciesMap;
    libhpx_cond_t dependenciesWaitCondition;
    libhpx_mutex_t dependenciesLock;
    floble_t dependenciesTimeNeuronWaitsFor;
  };
};

};  // algorithm

};  // neurox
