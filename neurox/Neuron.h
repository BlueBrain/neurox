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
    void sendSpikes(spike_time_t t); ///> fires AP, returns LCO for sent synapses
    size_t getSynapsesCount(); ///>number of netcons
    bool checkAPthresholdAndTransmissionFlag (floble_t v); ///> checks if AP threshold has been reached and whether spikes can be transmitted  (PreSynHelper->flag)

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
    void addSynapse(Synapse target);///> add hpx address of post-synaptic branch

    void informTimeDependantNeurons(floble_t t);
    void waitForTimeDependencyNeurons(floble_t t);
    void updateTimeDependencyTime(neuron_id_t srcGid, floble_t maxTimeAllowed);

  private:
    // the outgoing neurons:
    hpx_t synapsesMutex;   ///> mutex to protect variable 'synapses'
    std::vector<Synapse> synapses; ///> out-going synapse information
    bool synapsesTransmissionFlag; ///> PreSynHelper* psh -> flag (to flag whether spikes for a given AP have been sent or not

    ///the incoming neurons:
    std::map<neuron_id_t, floble_t> timeDependencies; ///> last known time of pre-synaptic neurons
    hpx_t timeDependenciesMutex;
    floble_t timeDependenciesMinTimeCached;
    void updateTimeDependenciesMinTimeChached();
}; //Neuron
}; //namespace neurox
