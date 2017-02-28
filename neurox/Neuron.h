#pragma once

#include "neurox/neurox.h"


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
        Synapse(hpx_t, floble_t, floble_t);
        hpx_t addr;        ///>address of destination
        floble_t nextNotificationTime; ///> next time this post-syn neuron needs to be informed of my actual time
        floble_t minDelay; ///>interval of notification in case of no spykes
                           ///(fastest Netcon from current neuron to dependant-neuron)
    };
    void sendSpikes(spike_time_t t); ///> fires AP, returns LCO for sent synapses
    void sendSteppingNotification(floble_t t); ///> inform my outgoing-connection neurons that I stepped
    void addSynapse(Synapse target);///> add hpx address of post-synaptic branch
    size_t getSynapseCount(); ///> get size of vector synapse

    //the incoming neuron connections:
    class TimeDependencies
    {
      public:
        TimeDependencies();
        ~TimeDependencies();

        //time dependency methods
        void waitForTimeDependencyNeurons(floble_t t);
        void updateTimeDependency(neuron_id_t srcGid, floble_t maxTimeAllowed, bool initialization = false);
        void updateDependenciesMinTimeCached();
      private:
        std::vector<std::pair<neuron_id_t, floble_t>> dependenciesVector; ///> pairs of <dependency-neuron-id, max-time-allowed> for pre-syn neurons
        std::map<neuron_id_t, floble_t*> dependenciesMap; ///> map to previous structure (pointing to vector values)
        floble_t dependenciesMinTimeCached; ///> last known minimal value of timeDependencies
        hpx_t dependenciesMutex; ///> protects dependenciesVector and dependenciesMinTimeCached

        hpx_t neuronWaitingLco;
        bool neuronWaitingFlag;
    } * timeDependencies;

  private:
    //the outgoing neuron connections:
    hpx_t synapsesMutex; ///> protects synapses
    std::vector<Synapse> synapses; ///> out-going synapse information
    bool synapsesTransmissionFlag; ///> PreSynHelper* psh -> flag (to flag whether spikes for a given AP have been sent or not

}; //Neuron
}; //namespace neurox
