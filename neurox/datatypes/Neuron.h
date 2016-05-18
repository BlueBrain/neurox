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

    //reporting vars
    int reportersCount;
    hpx_t * reporters;

    //outgoing synapses
    double APThreshold;
    int synapsesCount;      ///> number of outgoing synapses
    Synapse * synapses;     ///> outgoing Synapses

    static void registerHpxActions(); ///> Register all HPX actions
    static hpx_action_t initialize; ///> Initializes Neuron

  private:
    static int initialize_handler(const int gid, const hpx_t topBranch,
        double APThreshold, Synapse * synapses, size_t synapsesCount);
};
