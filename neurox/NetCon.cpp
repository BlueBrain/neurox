#include "neurox/neurox.h"
#include <cstring>

using namespace neurox;

NetConX::NetConX(){}

NetConX::NetConX(int mechType, offset_t mechInstance, floble_t delay,
                 floble_t * weights, unsigned short weightsCount, bool active)
    :mechType(mechType), weightsCount(weightsCount), mechInstance(mechInstance),
      delay(delay), active(active), weights(nullptr)
{
    if (weightsCount>0)
    {
        this->weights=new floble_t[weightsCount];
        memcpy(this->weights, weights, sizeof(floble_t)*weightsCount);
    }
}

void NetConX::deliver(floble_t tt, Branch* branch)
{
    if (!this->active) return
    getMechanismFromType(mechType)->callModFunction(branch, Mechanism::ModFunction::netReceive, this, tt);
}

NetConX::~NetConX()
{
    delete [] weights;
}
