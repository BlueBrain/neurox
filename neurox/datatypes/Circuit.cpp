#include "neurox/neurox.h"
#include "string.h"

Circuit::Circuit ():
  neuronsCount(0), multiSplit(0), neuronsAddr(HPX_NULL)
{}

Circuit::~Circuit()
{}


hpx_action_t Circuit::initialize = 0;
int Circuit::initialize_handler(const Circuit * circuit_new, const size_t size)
{
    //Make sure message arrived correctly, and pin memory
    hpx_t target = hpx_thread_current_target();
    uint64_t *local = NULL;
    if (!hpx_gas_try_pin(target, (void**) &local))
        return HPX_RESEND;

    //do the work
    circuit = new Circuit();
    memcpy(circuit, circuit_new, size);

    //unpin and return success
    hpx_gas_unpin(target);
    return HPX_SUCCESS;
}

void Circuit::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  initialize, initialize_handler, HPX_POINTER, HPX_SIZE_T);
}

