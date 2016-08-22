#pragma once

#include "neurox/Neurox.h"

namespace NeuroX
{

/**
 * @brief The Neuron class
 * Represents a neuron as its metadata and a tree-based orphology
 */
class Neuron
{
  public:
    Neuron() = delete;
    ~Neuron();

    //neuron metadata
    int id;					///> neuron global id
    hpx_t soma;		///> hpx address of the top compartment (soma)
    //TODO hpx_call_sync to the soma is 100x more expensive than function call
    //so I have to use the real soma or pointer instead (ask Luke)

    double t; ///> current time
    double dt; ///> time step size

    //outgoing synapses
    double APthreshold;  ///> Action Potential threshold
    std::vector<hpx_t> synapses;   ///> hpx address of post-synaptic recipient of synapse (neuron or branch)
    std::deque<hpx_t> synapsesLCO; ///> LCO for every AP sent

    void setupTreeMatrixMinimal(); ///>set_tree_matrix_minimal
    void callModFunction(Mechanism::ModFunction functionId); ///> Calls a MOD function on all mechanisms
    double getSomaVoltage(); ///> gets Soma voltage
    void fireActionPotential(); ///> fires AP, returns LCO for sent synapses
    void waitForSynapsesDelivery(int commStepSize); ///> waits for delivery of synapses

    static void registerHpxActions(); ///> Register all HPX actions
    static hpx_action_t init;         ///> Initializes Neuron
    static hpx_action_t addSynapseTarget;  ///> Inserts outgoing synapses (targets) in this Neuron
    static hpx_action_t broadcastNetCons;  ///> all branches broadcast to their pre-neuron ID their hpx addresses

  private:
    hpx_t synapsesMutex;   ///> mutex to protect variable 'synapses'

    static int init_handler(const int nargs, const void *args[], const size_t sizes[]); ///> HPX constructor
    static int addSynapseTarget_handler (const hpx_t * synapseTarget, const size_t size); ///>adds an outgoing Synapses
    static int broadcastNetCons_handler();
};

}
