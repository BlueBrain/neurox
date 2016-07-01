#pragma once

#include "neurox/Neurox.h"

namespace Neurox
{

/**
 * @brief The Mechanisms class
 * Stores unique mechanisms information and dependencies
 */
class Synapse
{
  public:
    Synapse();
    Synapse(const double weight, const double delay, const int mechType, const int mechInstance);
    ~Synapse();

    hpx_t target;  ///> address of target branch/neuron
    double weight; ///> synaptic weight
    double delay; ///synaptic delay

    bool operator<(const Synapse& rhs) const; ///> less-than operator

  private:
};

class SynapseOut
{
    hpx_t postNeuronAddr;
    int postNeuronId;
    double delay;
};

class SynapseIn
{
    int preNeuronId;
    int mechType;
    int mechInstance;
    double weight;
};

}
