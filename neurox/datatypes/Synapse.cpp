#include "neurox/Neurox.h"

using namespace Neurox;

NetConX::NetConX(){};

NetConX::NetConX(short int mechType, int mechInstance, double delay, double weight)
    :mechType(mechType), mechInstance(mechInstance), delay(delay), weight(weight) {};

Spike::Spike(const double deliveryTime, NetConX * synapse)
    :deliveryTime(deliveryTime), netcon(synapse) {};

Spike::~Spike(){};

bool Spike::operator<(const Spike& rhs) const
{
    return deliveryTime < rhs.deliveryTime;
}
