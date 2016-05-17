#include "neurox/neurox.h"
#include "string.h"

Brain brain; //global variable (defined in neurox.h)

Brain::Brain ():
  neuronsCount(0), multiSplit(0), neuronsAddr(HPX_NULL)
{}

hpx_action_t Brain::clear = 0;
int Brain::clear_handler()
{
    //Make sure message arrived correctly, and pin memory
    hpx_t target = hpx_thread_current_target();
    uint64_t *local = NULL;
    if (!hpx_gas_try_pin(target, (void**) &local))
        return HPX_RESEND;

    //do the work
    //for (int n=0; n<neuronsCount; n++)
    {
    //TODO delete everything
    }


    //unpin and return success
    hpx_gas_unpin(target);
    return HPX_SUCCESS;
}

hpx_action_t Brain::initialize = 0;
int Brain::initialize_handler(const Brain * brain_new, const size_t size)
{
    //Make sure message arrived correctly, and pin memory
    hpx_t target = hpx_thread_current_target();
    uint64_t *local = NULL;
    if (!hpx_gas_try_pin(target, (void**) &local))
        return HPX_RESEND;

    //do the work
    memcpy(&brain, brain_new, size);

    //unpin and return success
    hpx_gas_unpin(target);
    return HPX_SUCCESS;
}

void Brain::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  initialize, initialize_handler, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  clear, clear_handler);
}

