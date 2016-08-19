#include "neurox/Neurox.h"
#include <cstring>

using namespace NeuroX;

NetConX::NetConX(int mechType, int mechInstance, double delay, double * args, short int argsCount, bool active)
    :mechType(mechType), argsCount(argsCount), mechInstance(mechInstance), delay(delay), active(active), args(nullptr)
{
    if (argsCount>0)
    {
        this->args=new double[argsCount];
        memcpy(this->args, args, sizeof(double)*argsCount);
    }
}

void NetConX::deliver(double t, void* branch)
{
    if (!this->active) return
    getMechanismFromType(mechType)->callNetReceiveFunction(branch, this, t, 0);
}

NetConX::~NetConX()
{
    delete [] args;
}
