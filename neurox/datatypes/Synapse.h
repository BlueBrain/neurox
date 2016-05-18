#pragma once

#include "neurox/neurox.h"

/**
 * @brief The Mechanisms class
 * Stores unique mechanisms information and dependencies
 */
class Synapse
{
  public:
    Synapse(){};
    Synapse(const double weight, const double delay, const hpx_t target);
    ~Synapse();

    double weight; ///> synaptic weight
    double delay;  ///> delivery delay
    hpx_t target;  ///> hpx address of mechanism handling synapse on post-synaptic side
  private:
};
