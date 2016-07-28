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
    memcpy(Neurox::inputParams, data, sizeof(Neurox::Input::InputParams));
    neurox_hpx_unpin;
}

hpx_action_t setMechanisms = 0;
int setMechanisms_handler(const int nargs, const void *args[], const size_t sizes[])
{
    /**
     * nargs=3 where:
     * args[0] = array of all mechanisms info
     * args[1] = array of all mechanisms dependencies
     * args[2] = array of all mechanisms names (sym)
     */

    neurox_hpx_pin(uint64_t);
    assert(nargs==3);

    if (mechanisms!=nullptr)
        delete [] Neurox::mechanisms;

    mechanismsCount = sizes[0]/sizeof(Mechanism);
    mechanisms = new Neurox::Mechanism[mechanismsCount];

    int offsetDependencies=0;
    int offsetSym=0;
    for (int m=0; m<mechanismsCount; m++)
    {
        Mechanism & mech = ((Mechanism*) args[0])[m];
        int * dependencies = &((int*) args[1])[offsetDependencies];
        char * sym = &((char*) args[2])[offsetSym];

        mechanisms[m] = Mechanism (mech.type, mech.dataSize, mech.pdataSize,
                                   mech.isArtificial, mech.pntMap, mech.isIon,
                                   mech.conci, mech.conco, mech.charge,
                                   mech.symLength, mech.sym,
                                   mech.dependenciesCount, dependencies);
        offsetDependencies +=  mech.dependenciesCount;
        offsetSym += mech.symLength;
    }
    neurox_hpx_unpin;
}

void registerHpxActions()
{
    neurox_hpx_register_action(1,Neurox::setInputParams);
    neurox_hpx_register_action(2,Neurox::setMechanisms);
}

};
