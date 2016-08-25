#pragma once

#include "neurox/neurox.h"

namespace neurox
{

/**
 * @brief The Neuron class
 * Represents a neuron as its metadata and a tree-based orphology
 */
class Neuron
{
  public:
    Neuron() = delete;
    Neuron(int neuronId, double APthreshold);
    ~Neuron();

    int id;              ///> neuron global id
    double APthreshold;  ///> Action Potential threshold
    void setupTreeMatrixMinimal(); ///>set_tree_matrix_minimal
    void fireActionPotential(double t); ///> fires AP, returns LCO for sent synapses
    void waitForSynapsesDelivery(int commStepSize); ///> waits for delivery of synapses
    size_t getNetConsCount(); ///>number of netcons
    void addSynapseTarget(hpx_t target);///> add hpx address of post-synaptic branch

  private:
    hpx_t synapsesMutex;   ///> mutex to protect variable 'synapses'
    std::vector<hpx_t> synapses;   ///> hpx address of post-synaptic recipient of synapse (as a hpx address of a branch)
    std::deque<hpx_t> synapsesLCO; ///> LCO for every AP sent
};
}
