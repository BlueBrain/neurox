#pragma once

#include "neurox/neurox.h"

/**
 * @brief The Mechanisms class
 * Stores unique mechanisms information and dependencies
 */
class Synapse
{
  public:
    Synapse()=delete;
    Synapse(const double weight, const double delay, const int mechType, const int mechInstance);
    ~Synapse();

    double weight; ///> synaptic weight
    double delay;  ///> delivery delay
    int mechType; ///> mechanism type
    int mechInstance; ///> instance of this mechanism, according to the mech type
  private:
};
