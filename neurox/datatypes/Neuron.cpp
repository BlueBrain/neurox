#include "neurox/neurox.h"
#include "string.h"

Neuron::~Neuron()
{
    delete [] synapses;
}

hpx_action_t Neuron::initialize = 0;
int Neuron::initialize_handler(const int gid, const hpx_t topBranch,
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
    neuron->synapses = new Symbol[synapsesCount];
    memcpy(neuron->synapses, synapses, sizeof(Synapse)*synapsesCount);

    //unpin and return success
    hpx_gas_unpin(neuron_addr);
    return HPX_SUCCESS;
}

void Neuron::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  initialize, initialize_handler, HPX_INT, HPX_ADDR, HPX_DOUBLE, HPX_POINTER, HPX_SIZE_T);
}


