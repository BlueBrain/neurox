#include "neurox/neurox.h"
#include <cstring>

Brain brain; //global variable (defined in neurox.h)

Brain::Brain ():
  neuronsCount(0), neuronsAddr(HPX_NULL)
{

}

Brain::~Brain()
{
    delete [] mechanisms;
}

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
int Brain::initialize_handler(const Brain * brain_new,
           const Mechanism * mechanisms, const int mechanismsCount,
           const int * mechDependencies, const int totalDependenciesCount)
{
    //Make sure message arrived correctly, and pin memory
    hpx_t target = hpx_thread_current_target();
    uint64_t *local = NULL;
    if (!hpx_gas_try_pin(target, (void**) &local))
        return HPX_RESEND;

    //blind copy of brain and mechanisms info
    memcpy(&brain, &brain_new, sizeof(Brain));
    memcpy(&brain.mechanisms, &brain_new->mechanisms, mechanismsCount*sizeof(Mechanism));

    //add mechanisms information
    int acc=0;
    for (int m=0; m<mechanismsCount; m++)
    {
        int dependenciesCount = brain.mechanisms[m].dependenciesCount;
        brain.mechanisms[m].dependencies = new int[dependenciesCount];
        std::memcpy(brain.mechanisms[m].dependencies, &mechDependencies[acc], dependenciesCount*sizeof(int));
        acc+=dependenciesCount;
    }

    //clean up
    delete [] mechanisms;
    delete [] mechDependencies;
    delete brain_new;

    //unpin and return success
    hpx_gas_unpin(target);
    return HPX_SUCCESS;
}

void Brain::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  initialize, initialize_handler, HPX_POINTER,  HPX_POINTER, HPX_INT, HPX_POINTER, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  clear, clear_handler);
}

