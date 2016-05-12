#include "neurox/neurox.h"
#include "string.h"

//global variable
Brain brain;

Brain::Brain ():
  neuronsCount(0), multiSplit(0), neuronsAddr(HPX_NULL)
{}

Brain::~Brain()
{
    //TODO
}

hpx_action_t Brain::initialize = 0;
int Brain::initialize_handler(const Brain * circuit_new, const size_t size)
{
    //Make sure message arrived correctly, and pin memory
    hpx_t target = hpx_thread_current_target();
    uint64_t *local = NULL;
    if (!hpx_gas_try_pin(target, (void**) &local))
        return HPX_RESEND;

    //do the work
    memcpy(&circuit, circuit_new, size);

    //unpin and return success
    hpx_gas_unpin(target);
    return HPX_SUCCESS;
}

void Brain::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  initialize, initialize_handler, HPX_POINTER, HPX_SIZE_T);
}

