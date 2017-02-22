#include "neurox/neurox.h"
#include <cstring>

using namespace neurox;

NetConX::NetConX(){}

NetConX::NetConX(int mechType, offset_t mechInstance, floble_t delay,
                 offset_t weightIndex, unsigned short weightsCount, bool active)
    :mechType(mechType), weightsCount(weightsCount), mechInstance(mechInstance),
      delay(delay), weightIndex(weightIndex), active(active)
{}

NetConX::~NetConX()
{}

void NetConX::deliver(floble_t tt, Branch* branch)
{
    if (!this->active) return
    getMechanismFromType(mechType)->callModFunction(branch, Mechanism::ModFunction::netReceive, this, tt);
}

