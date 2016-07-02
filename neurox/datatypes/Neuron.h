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
    double APthreshold;      ///> Action Potential threshold
    int synapsesCount;       ///> number of outgoing synapses
    SynapseOut * synapses;   ///> outgoing synapses

    static void registerHpxActions(); ///> Register all HPX actions
    static void setupTreeMatrixMinimal(Neuron * local); ///>set_tree_matrix_minimal
    static hpx_t fireActionPotential(Neuron * local); ///> fires AP, returns LCO for sent synapses
    static hpx_action_t finitialize;  ///> finitialize.c
    static hpx_action_t init;         ///> Initializes Neuron
    static hpx_action_t addOutgoingSynapses;  ///> Inserts outgoing synapses (targets) in this Neuron


  private:
    static int finitialize_handler(); ///> initialize.c
    static int init_handler(const int gid, const hpx_t topBranch, double APthreshold); ///> HPX constructor
    static int addOutgoingSynapses_handler (const SynapseOut * synapsesOut, const size_t size); ///>Sets all outgoing Synapses
};

}
