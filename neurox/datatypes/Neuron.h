#pragma once

#include "neurox/neurox.h"

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
    hpx_t topBranch;		///> hpx address of the top compartment (soma)
    
    //from NrnThread
    double cj; ///<1st or 2nd order solver ... (?)
    
    //reporting vars
    int reportersCount;
    hpx_t * reporters;

    //outgoing synapses
    double APThreshold;
    int synapsesCount;      ///> number of outgoing synapses
    Synapse * synapses;     ///> outgoing Synapses

    static void registerHpxActions(); ///> Register all HPX actions
    static hpx_action_t init;         ///> Initializes Neuron
    static hpx_action_t finitialize;  ///> finitialize.c::nrn_finitialize

  private:
    static int init_handler(const int gid, const hpx_t topBranch,
        double APThreshold, Synapse * synapses, size_t synapsesCount);

    static int finitialize_handler(); ///finitialize.c
};
