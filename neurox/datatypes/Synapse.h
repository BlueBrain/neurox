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

    double weight; ///> synaptic weight
    double delay; ///synaptic delay
    int mechType; ///> mechanism type
    int mechInstance; ///> instance of this mechanism, according to the mech type
    double deliveryTime;  ///> last delivery time (time of AP + delay)

    bool operator<(const Synapse& rhs) const; ///> less-than operator

  private:
};

}
