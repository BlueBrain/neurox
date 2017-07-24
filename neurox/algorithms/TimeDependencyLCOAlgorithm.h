#pragma once
#include "neurox.h"

#define DERIVED_CLASS_NAME TimeDependencyLCOAlgorithm

using namespace neurox;

namespace neurox
{

namespace algorithms
{

class DERIVED_CLASS_NAME : public Algorithm
{
  public:
    DERIVED_CLASS_NAME();
    ~DERIVED_CLASS_NAME();

    const AlgorithmType getType() override;
    const char* getTypeString() override;

    void Init() override;
    void Clear() override;
    double Launch() override;

    void StepBegin(Branch*) override;
    void StepEnd(Branch*, hpx_t) override;
    void Run(Branch*, const void*) override;
    hpx_t SendSpikes(Neuron* b, double tt, double t) override;
    void AfterReceiveSpikes( Branch *local, hpx_t target, neuron_id_t preNeuronId,
                             spike_time_t spikeTime, spike_time_t maxTime) override;


    ///controls time-dependencies from incoming neuron connections
    class TimeDependencies : public AlgorithmMetaData
    {
      public:
        TimeDependencies();
        ~TimeDependencies();

        void WaitForTimeDependencyNeurons(floble_t t, floble_t dt, int gid);
        void SendSteppingNotification(floble_t t, floble_t dt, int gid, std::vector<Neuron::Synapse*> & synapses); ///> inform my outgoing-connection neurons that I stepped
        void UpdateTimeDependency(neuron_id_t srcGid, floble_t dependencyNotificationTime, neuron_id_t myGid = -1, bool initialization = false);
        floble_t GetDependenciesMinTime();
        size_t GetDependenciesCount();
        void IncreseDependenciesTime(floble_t t);

        static floble_t notificationIntervalRatio; ///> ratio of notification interval (0,1]
        static double teps; ///>time-epsilon to correct wrong delivery of events due to floating point rounding

      private:
        std::map<neuron_id_t, floble_t> dependenciesMap; ///> map to previous structure (pointing to vector values)
        libhpx_cond_t dependenciesWaitCondition;
        libhpx_mutex_t dependenciesLock;
        floble_t dependenciesTimeNeuronWaitsFor;
    };
};

}; //algorithm

}; //neurox
