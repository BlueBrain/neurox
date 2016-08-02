#include "neurox/Neurox.h"
#include <cstring>

using namespace Neurox;

NetConX::NetConX(){};

NetConX::NetConX(int mechType, int mechInstance, double delay, double * args, short int argsCount, char active)
    :mechType(mechType), argsCount(argsCount), mechInstance(mechInstance), delay(delay), active(active), args(nullptr)
{
    if (argsCount>0)
    {
        this->args=new double[argsCount];
        memcpy(this->args, args, sizeof(double)*argsCount);
    }
}

NetConX::~NetConX()
{
    delete [] args;
}
