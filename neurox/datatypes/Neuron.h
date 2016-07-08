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
    static hpx_t fireActionPotential(Neuron * local); ///> fires AP, returns LCO for sent synapses
    static hpx_action_t finitialize;  ///> finitialize.c
    static hpx_action_t init;         ///> Initializes Neuron
    static hpx_action_t addSynapseTarget;  ///> Inserts outgoing synapses (targets) in this Neuron

  private:
    hpx_t synapsesMutex;   ///> mutex to protect the memory access to synapses vector

    static int finitialize_handler(); ///> initialize.c
    static int init_handler(const int gid, const hpx_t topBranch, double APthreshold); ///> HPX constructor
    static int addSynapseTarget_handler (const hpx_t synapseTarget); ///>adds an outgoing Synapses
};

}
