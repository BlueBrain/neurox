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
    Synapse(const double weight, const double delay, const hpx_t target, const int mechOffset, const int mechInstance);
    ~Synapse();

    double weight; ///> synaptic weight
    double delay;  ///> delivery delay
    hpx_t target;  ///> hpx address of branch holding target mech
    int mechOffset;///> offset of mechanisms in this branch
    int mechInstance; ///> instance of this mechanism, according to the mech type
  private:
};
