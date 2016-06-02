#include "neurox/Neurox.h"

using namespace Neurox;

Synapse::Synapse(const double weight, const double delay, const int mechType, const int mechInstance)
    :weight(weight), delay(delay), mechType(mechType), mechInstance(mechInstance), deliveryTime(0) {};

Synapse::~Synapse(){};

bool Synapse::operator<(const Synapse& rhs) const
{
    return deliveryTime < rhs.deliveryTime;
}
