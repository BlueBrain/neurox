#include "neurox/Neurox.h"

using namespace Neurox;

SynapseOut::SynapseOut(){};

SynapseOut::SynapseOut(hpx_t postNeuronAddr, int postNeuronId, double delay)
    :postNeuronAddr(postNeuronAddr), postNeuronId(postNeuronId), delay(delay){};

SynapseIn::SynapseIn(){};

SynapseIn::SynapseIn(int preNeuronId, int mechType, int mechInstance, double weight)
    :preNeuronId(preNeuronId), mechType(mechType), mechInstance(mechInstance), weight(weight) {};

Spike::Spike(const double deliveryTime, SynapseIn * synapse)
    :deliveryTime(deliveryTime), synapse(synapse) {};

Spike::~Spike(){};

bool Spike::operator<(const Spike& rhs) const
{
    return deliveryTime < rhs.deliveryTime;
}
