#include "neurox/neurox.h"
#include <cstring>

using namespace neurox;

NetConX::NetConX(){}

NetConX::NetConX(int mechType, offset_t mechInstance, floble_t delay,
                 floble_t * args, unsigned short argsCount, bool active)
    :mechType(mechType), argsCount(argsCount), mechInstance(mechInstance),
      delay(delay), active(active), args(nullptr)
{
    if (argsCount>0)
    {
        this->args=new floble_t[argsCount];
        memcpy(this->args, args, sizeof(floble_t)*argsCount);
    }
}

void NetConX::deliver(floble_t t, Branch* branch)
{
    if (!this->active) return
    getMechanismFromType(mechType)->callNetReceiveFunction(branch, this, t, 0);
}

NetConX::~NetConX()
{
    delete [] args;
}
