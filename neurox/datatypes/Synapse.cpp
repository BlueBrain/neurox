#include "neurox/Neurox.h"

using namespace Neurox;

Synapse::Synapse(const hpx_t target, const double weight, const double delay, const int mechType, const int mechInstance)
    :target(target), weight(weight), delay(delay) {};

Synapse::Synapse(){};

Synapse::~Synapse(){};

bool Synapse::operator<(const Synapse& rhs) const
{
    return deliveryTime < rhs.deliveryTime;
}
