#pragma once

#include "neurox/Neurox.h"

namespace NeuroX
{

/**
 * @brief The Spike delivery class
 * Equivalent to NRNMPI_Spike in Coreneuron (nrnmpi.h)
 * Stores the information of a spike to be delivered (delivery time and Synaptic input info)
 */
class Spike
{
  public:
    Spike() = delete;
    Spike(const double deliveryTime, NetConX * netcon);
    ~Spike();

    double deliveryTime;  ///> delivery time of synapse
    NetConX * netcon;     ///> synapse to be delivered
    bool operator<(const Spike& rhs) const; ///> less-than operator

  private:
};

}
