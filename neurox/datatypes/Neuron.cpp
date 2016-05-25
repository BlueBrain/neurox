#include "neurox/neurox.h"
#include "string.h"

Neuron::~Neuron()
{
    delete [] synapses;
}

hpx_action_t Neuron::init = 0;
int Neuron::init_handler(const int gid, const hpx_t topBranch,
    double APThreshold, Synapse * synapses, size_t synapsesCount)
{
    //Make sure message arrived correctly, and pin memory
    hpx_t neuron_addr = hpx_thread_current_target();
    Neuron * neuron = NULL;
    if (!hpx_gas_try_pin(neuron_addr, (void**) &neuron))
        return HPX_RESEND;

    //copy information
    neuron->topBranch=topBranch;
    neuron->id=gid;
    neuron->synapsesCount=synapsesCount;
    neuron->APThreshold=APThreshold;
    neuron->synapses = new Synapse[synapsesCount];
    memcpy(neuron->synapses, synapses, sizeof(Synapse)*synapsesCount);

    //unpin and return success
    hpx_gas_unpin(neuron_addr);
    return HPX_SUCCESS;
}

hpx_action_t Neuron::finitialize = 0;
int Neuron::finitialize_handler()
{
    //Make sure message arrived correctly, and pin memory
    hpx_t neuron_addr = hpx_thread_current_target();
    Neuron * neuron = NULL;
    if (!hpx_gas_try_pin(neuron_addr, (void**) &neuron))
        return HPX_RESEND;

    //set up by finitialize.c:nrn_finitialize() -> fadvance_core.c:dt2thread()
    neuron->cj = inputParams->secondorder ? 2.0/inputParams->dt : 1.0/inputParams->dt;

    //call it for top branch
    hpx_call_sync(neuron->topBranch, Branch::finitialize);}

    //unpin and return success
    hpx_gas_unpin(neuron_addr);
    return HPX_SUCCESS;
}

void Neuron::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  init, init_handler, HPX_INT, HPX_ADDR, HPX_DOUBLE, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  finitialize, finitialize_handler);
}


