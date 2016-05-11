#include "neurox/neurox.h"
#include "string.h"

hpx_action_t Neuron::initialize = 0;
int Neuron::initialize_handler(const Neuron * inputNeuron,  const size_t size)
{
    //Make sure message arrived correctly, and pin memory
    hpx_t neuron_addr = hpx_thread_current_target();
    Neuron * neuron = NULL;
    if (!hpx_gas_try_pin(neuron_addr, (void**) &neuron))
        return HPX_RESEND;

    //copy header information
    memcpy(neuron, inputNeuron, size);

    //unpin and return success
    hpx_gas_unpin(neuron_addr);
    return HPX_SUCCESS;
}


void Neuron::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  initialize, initialize_handler, HPX_POINTER, HPX_SIZE_T);
}


