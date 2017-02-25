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
    size_t getSynapsesCount(); ///>number of netcons
    void addSynapseTarget(hpx_t target);///> add hpx address of post-synaptic branch
    bool checkAPthresholdAndTransmissionFlag (floble_t v); ///> checks if AP threshold has been reached and whether spikes can be transmitted  (PreSynHelper->flag)

    class SlidingTimeWindow
    {
      public:
        SlidingTimeWindow()=delete;
        SlidingTimeWindow(floble_t windowSize);
        ~SlidingTimeWindow();

        hpx_t dependant_lco; ///> lco of my dependant
        hpx_t dependencies_lco; ///> lco of my dependencies
        floble_t dependencies_t; ///> last-known time of my dependencies
        floble_t windowSize;

        void informTimeDependants(floble_t t);
        void waitForTimeDependencies(floble_t t);

        static hpx_action_t initDependencies; ///> Initializes dependencies
        static hpx_action_t setDependant;

        static int initDependencies_handler();
        static int setDependant_handler(const hpx_t * dependant_ptr, const size_t);
    } * slidingTimeWindow;

    static void registerHpxActions();

  private:
    hpx_t synapsesMutex;   ///> mutex to protect variable 'synapses'
    std::vector<hpx_t> synapsesTargets;   ///> hpx address of post-synaptic synapse targets
    bool synapsesTransmissionFlag; ///> PreSynHelper* psh -> flag (to flag whether spikes for a given AP have been sent or not
};
}
