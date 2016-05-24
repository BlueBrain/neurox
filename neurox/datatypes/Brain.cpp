#include "neurox/neurox.h"
#include <cstring>

Brain * brain = nullptr; //global variable (defined in neurox.h)

Brain::~Brain()
{
    delete [] mechsTypes;
}

Brain::Brain(const int neuronsCount,
             const hpx_t neuronsAddr, const Mechanism * mechanisms,
             const size_t mechanismsCount, const int * mechDependencies)
    : neuronsCount(neuronsCount), neuronsAddr(neuronsAddr), mechsTypesCount(mechanismsCount)
{
    //add mechanisms information
    int offset=0;
    this->mechsTypes = new Mechanism[mechanismsCount];
    for (int m=0; m<mechanismsCount; m++)
    {
        this->mechsTypes[m]=Mechanism(mechanisms[m].dataSize, mechanisms[m].pdataSize,
                                      mechanisms[m].dependenciesCount, mechanisms[m].pntMap,
                                      mechanisms[m].isArtificial, &mechDependencies[offset]);
        offset += mechanisms[m].dependenciesCount;
    }
};

hpx_action_t Brain::clear = 0;
int Brain::clear_handler()
{
    //Make sure message arrived correctly, and pin memory
    hpx_t target = hpx_thread_current_target();
    uint64_t *local = NULL;
    if (!hpx_gas_try_pin(target, (void**) &local))
        return HPX_RESEND;

    delete brain;

    //unpin and return success
    hpx_gas_unpin(target);
    return HPX_SUCCESS;
}

hpx_action_t Brain::initialize = 0;
int Brain::initialize_handler(const int neuronsCount,
           const hpx_t neuronsAddr, const Mechanism * mechanisms,
           const size_t mechanismsCount, const int * mechDependencies)
{
    //Make sure message arrived correctly, and pin memory
    hpx_t target = hpx_thread_current_target();
    uint64_t *local = NULL;
    if (!hpx_gas_try_pin(target, (void**) &local))
        return HPX_RESEND;

    //initialize global variable brain
    brain = new Brain(neuronsCount, neuronsAddr, mechanisms, mechanismsCount, mechDependencies);

    //clean up
    delete [] mechanisms;
    delete [] mechDependencies;

    //unpin and return success
    hpx_gas_unpin(target);
    return HPX_SUCCESS;
}

void Brain::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  initialize, initialize_handler, HPX_INT, HPX_ADDR, HPX_POINTER, HPX_SIZE_T, HPX_POINTER);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  clear, clear_handler);
}

