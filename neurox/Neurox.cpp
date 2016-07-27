#include "neurox/Neurox.h"
#include <cstring>

namespace Neurox
{

int neuronsCount=-1;
hpx_t neuronsAddr = HPX_NULL;
int mechanismsCount=-1;
Mechanism * mechanisms = nullptr;
Input::InputParams * inputParams = nullptr;

hpx_action_t setInputParams = 0;
int setInputParams_handler(const Input::InputParams * data, const size_t)
{
    neurox_hpx_pin(uint64_t);
    if (mechanisms!=nullptr)
        delete [] Neurox::inputParams;

    inputParams = new Neurox::Input::InputParams();
    memcpy(Neurox::inputParams, data, sizeof(inputParams));
    neurox_hpx_unpin;
}

hpx_action_t setMechanisms = 0;
int setMechanisms_handler(const int nargs, const void *args[], const size_t[])
{
    /**
     * nargs=4 where:
     * args[0] = number of mechanisms
     * args[1] = array of all mechanisms info
     * args[2] = number of dependencies per mechanism
     * args[3] = array of all mechanisms dependencies
     */

    neurox_hpx_pin(uint64_t);
    assert(nargs==4);

    if (mechanisms!=nullptr)
        delete [] Neurox::mechanisms;

    //Copy basic mechanisms data
    mechanismsCount = *((int*) args[0]);
    mechanisms = new Neurox::Mechanism[mechanismsCount];
    memcpy(mechanisms, args[1], mechanismsCount*sizeof(Mechanism));

    //Copy dependencies data
    int offset=0;
    for (int m=0; m<mechanismsCount; m++)
    {
        int dependenciesCount = ((int*) args[2])[m];
        int * dependencies = &((int*) args[3])[offset];

        mechanisms[m].dependenciesCount = dependenciesCount;
        mechanisms[m].dependencies = new int[dependenciesCount];
        memcpy(mechanisms[m].dependencies, dependencies, dependenciesCount*sizeof(int));
        offset += dependenciesCount;
    }
    neurox_hpx_unpin;
}

void registerHpxActions()
{
    neurox_hpx_register_action(1,Neurox::setInputParams);
    neurox_hpx_register_action(3,Neurox::setMechanisms);
}

};
