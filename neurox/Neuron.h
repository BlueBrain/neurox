#pragma once

#include "neurox/neurox.h"
#include "libhpx/libhpx.h"
#include <deque>

namespace neurox
{

/**
 * @brief The Neuron class
 * Represents the soma structure and outgoing spikes network of a neuron
 */
class Neuron
{
  public:
    Neuron() = delete;
    ~Neuron();

    Neuron(neuron_id_t neuronId, floble_t APthreshold, int thvar_index);

    static const double teps; ///>time-epsilon to correct wrong delivery of events due to floating point rounding

    neuron_id_t gid;     ///> neuron global id
    floble_t threshold;  ///> Action Potential threshold (PreSyn->_threshold)
    int thvar_index;     ///> index in voltages array og value to be compared against AP threshold (PreSyn->thvar_index_)
    floble_t refractoryPeriod; ///> refractory period

    void setupTreeMatrixMinimal(); ///>set_tree_matrix_minimal
    bool checkAPthresholdAndTransmissionFlag (floble_t v); ///> checks if AP threshold has been reached and whether spikes can be transmitted  (PreSynHelper->flag)

    // the outgoing neurons:
    class Synapse
    {
      public:
        Synapse()=delete;
        Synapse(hpx_t addr, floble_t minDelay, int destinationGid=-1);
        ~Synapse();
        hpx_t addr;        ///>address of destination
#ifndef NDEBUG
        int destinationGid;
#endif
        floble_t nextNotificationTime; ///> next time this post-syn neuron needs to be informed of my actual time
        floble_t minDelay; ///>interval of notification in case of no spykes
                           ///(fastest Netcon from current neuron to dependant-neuron)
        hpx_t previousSpikeLco; ///>lco controlling spikes delivery
    };

    hpx_t commStepAllSpikesLco; ///> LCO for all spikes of previous Comm Step (for fixed step methods and debug)
    static const int commStepSize; ///> Fixed communication step size

    void sendSpikes(floble_t t, floble_t dt); ///> fires AP, returns LCO for sent synapses
    void addSynapse(Synapse * target);///> add hpx address of post-synaptic branch
    size_t getSynapseCount(); ///> get size of vector synapse

    //the incoming neuron connections:
    class TimeDependencies
    {
      public:
        TimeDependencies();
        ~TimeDependencies();

        //time dependency methods
        void waitForTimeDependencyNeurons(floble_t t, floble_t dt, int gid);
        void sendSteppingNotification(floble_t t, floble_t dt, int gid, std::vector<Synapse*> & synapses); ///> inform my outgoing-connection neurons that I stepped
        void updateTimeDependency(neuron_id_t srcGid, neuron_id_t myGid, floble_t dependencyNotificationTime, bool initialization = false);
        floble_t getDependenciesMinTime();
        size_t getDependenciesCount();

        ///> how often to communicate stepping notification (1: every 'minDelay', 0.5 every half 'minDelay)
        static const floble_t notificationIntervalRatio;

      private:
        std::map<neuron_id_t, floble_t> dependenciesMap; ///> map to previous structure (pointing to vector values)
        libhpx_cond_t dependenciesWaitCondition;
        libhpx_mutex_t dependenciesLock;
        floble_t dependenciesTimeNeuronWaitsFor;

    } * timeDependencies;

    std::vector<Synapse*> synapses; ///> out-going synapse information
  private:
    //the outgoing neuron connections:
    hpx_t synapsesMutex; ///> protects synapses
    bool synapsesTransmissionFlag; ///> PreSynHelper* psh -> flag (to flag whether spikes for a given AP have been sent or not

}; //Neuron
}; //namespace neurox
