#include "neurox/Neurox.h"

using namespace Neurox;

namespace Neurox
{
    int neuronsCount=-1; 
    hpx_t neuronsAddr=HPX_NULL;
    int mechanismsCount=-1;
    Mechanism * mechanisms=NULL;
    Input::InputParams * inputParams = nullptr;

    hpx_action_t setInputParams = 0;
    int setInputParams_handler(const Input::InputParams * inputParams, const size_t size)
    {
        neurox_hpx_pin(uint64_t);
        if (Neurox::inputParams==NULL)
            delete Neurox::inputParams;

        Neurox::inputParams = new Input::InputParams();
        memcpy(Neurox::inputParams, inputParams, size);
        neurox_hpx_unpin;
    };

    hpx_action_t setMechanisms = 0;
    int setMechanisms_handler(const Mechanism * mechanisms, const size_t count)
    {
        neurox_hpx_pin(uint64_t);
        if (Neurox::mechanisms==NULL)
            delete [] Neurox::mechanisms;

        Neurox::mechanisms = new Mechanism[count];
        memcpy(Neurox::mechanisms, mechanisms, sizeof(Mechanism)*count);
        neurox_hpx_unpin;
    }
};
