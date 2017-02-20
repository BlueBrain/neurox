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
    void setupTreeMatrixMinimal(); ///>set_tree_matrix_minimal
    void sendSpikes(spike_time_t t); ///> fires AP, returns LCO for sent synapses
    void waitForSynapsesDelivery(int commStepSize); ///> waits for delivery of synapses
    size_t getSynapsesCount(); ///>number of netcons
    void addSynapseTarget(hpx_t target);///> add hpx address of post-synaptic branch
    bool checkAPthresholdAndTransmissionFlag (floble_t v); ///> checks if AP threshold has been reached and whether spikes can be transmitted  (PreSynHelper->flag)

  private:
    hpx_t synapsesMutex;   ///> mutex to protect variable 'synapses'
    std::vector<hpx_t> synapsesTargets;   ///> hpx address of post-synaptic synapse targets
    std::deque<hpx_t> synapsesLCO; ///> LCO for every AP sent
    bool synapsesTransmissionFlag; ///> PreSynHelper* psh -> flag (to flag whether spikes for a given AP have been sent or not
};
}
