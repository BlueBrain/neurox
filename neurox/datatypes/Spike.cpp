#include "neurox/Neurox.h"

using namespace NeuroX;

Spike::Spike(const double deliveryTime, NetConX * synapse)
    :deliveryTime(deliveryTime), netcon(synapse) {}

Spike::~Spike(){}

bool Spike::operator<(const Spike& rhs) const
{
    return deliveryTime < rhs.deliveryTime;
}
