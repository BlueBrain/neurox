#pragma once

#include "neurox/Neurox.h"

namespace Neurox
{

/**
 * @brief The SynapseOut class
 * Includes the synaptic information at the pre-synaptic neuron
 */
class SynapseOut
{
  public:
    SynapseOut();
    SynapseOut(hpx_t postNeuronAddr, int postNeuronId, double delay);

    hpx_t postNeuronAddr; ///> post-synaptic neuron address
    int postNeuronId;     ///> post-synaptic neuron id
    double delay;         ///> synaptic delay (soma-bouton distance + transmitters release delay)
};

/**
 * @brief The SynapseIn class
 * Includes the synaptic information at the post-synaptic neuron
 */
class SynapseIn
{
  public:
    SynapseIn();
    SynapseIn(int preNeuronId, int mechType, int mechInstance, double weight);

    int preNeuronId;      ///> pre-synaptic neuron id
    int mechType;         ///> mechanism type associated with this synapse
    int mechInstance;     ///> mechanism instance, from the mechanism type
    double weight;        ///> synaptic weight (0 means disabled)
};

/**
 * @brief The Spike delivery class
 * Stores the information of a spike to be delivered (delivery time and Synaptic input info)
 */
class Spike
{
  public:
    Spike() = delete;
    Spike(const double deliveryTime, SynapseIn * synapse);
    ~Spike();

    double deliveryTime;   ///> delivery time of synapse
    SynapseIn * synapse;   ///> synapse to be delivered
    bool operator<(const Spike& rhs) const; ///> less-than operator

  private:
};

}
