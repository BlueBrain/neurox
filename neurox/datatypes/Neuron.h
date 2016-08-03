#pragma once

#include "neurox/Neurox.h"

namespace Neurox
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

    double t; ///> current time
    double dt; ///> time step size

    //from NrnThread
    double cj; ///<1st or 2nd order solver ... (?)

    //outgoing synapses
    double APthreshold;  ///> Action Potential threshold
    std::vector<hpx_t> synapses;   ///> hpx address of post-synaptic recipient of synapse (neuron or branch)

    static void registerHpxActions(); ///> Register all HPX actions
    static void setupTreeMatrixMinimal(Neuron * local); ///>set_tree_matrix_minimal
    static void callMechanismsModFunction(Neuron * local, Mechanism::ModFunction functionId); ///> Calls MOD functions on all mechanisms
    static hpx_t fireActionPotential(Neuron * local); ///> fires AP, returns LCO for sent synapses
    static hpx_action_t finitialize;  ///> finitialize.c
    static hpx_action_t init;         ///> Initializes Neuron
    static hpx_action_t addSynapseTarget;  ///> Inserts outgoing synapses (targets) in this Neuron

  private:
    hpx_t synapsesMutex;   ///> mutex to protect the memory access to synapses vector

    static int init_handler(const int nargs, const void *args[], const size_t sizes[]); ///> HPX constructor
    static int addSynapseTarget_handler (const hpx_t * synapseTarget, const size_t size); ///>adds an outgoing Synapses
    static int finitialize_handler(); ///> initialize.c
};

}
