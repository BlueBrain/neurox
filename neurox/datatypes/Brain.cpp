#include "neurox/neurox.h"
#include <cstring>

Brain * brain = nullptr; //global variable (defined in neurox.h)

Brain::~Brain()
{
    delete [] mechsTypes;
}

Brain::Brain(const int neuronsCount,
             const hpx_t neuronsAddr, const Mechanism * mechanisms,
             const size_t mechsTypesCount, const int * mechDependencies)
    : neuronsCount(neuronsCount), neuronsAddr(neuronsAddr), mechsTypesCount(mechsTypesCount)
{
    //add mechanisms information
    int offset=0;
    this->mechsTypes = new Mechanism[mechsTypesCount];
    for (int m=0; m<mechsTypesCount; m++)
    {
        this->mechsTypes[m]=Mechanism(m, mechanisms[m].dataSize, mechanisms[m].pdataSize,
                                      mechanisms[m].dependenciesCount, mechanisms[m].pntMap,
                                      mechanisms[m].isArtificial, &mechDependencies[offset],
                                      mechanisms[m].isIon, mechanisms[m].conci, mechanisms[m].conco, mechanisms[m].charge //for ions
                                      );
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

hpx_action_t Brain::finitialize = 0;
int Brain::finitialize_handler()
{
    //Make sure message arrived correctly, and pin memory
    hpx_t brain_addr = hpx_thread_current_target();
    Brain *brain = NULL;
    if (!hpx_gas_try_pin(brain_addr, (Brain**) &brain))
        return HPX_RESEND;

    hpx_par_for_sync(
       [&] (int i, void*) { hpx_call_sync(brain->getNeuronAddr(i), Neuron::finitialize);},
       0, brain->neuronsCount, NULL);

    //TODO detect min delay (must be divisable by dt)
    //in netpar.c: usable_mindelay_ = floor(mindelay_ * dt1_ + 1e-9) * dt;
    //probably should be solver-depentend eg Solver::init(...)

    //TODO nrn_deliver_events()

    //BEFORE_INITIAL
    //for (int m=0; m<brain->mechsTypesCount; m++)
        //beforeAfterFunctions[BEFORE_INITIAL](NrnThread nt, Memb_list * ml, type);

    //INITIAL
    //for (int m=0; m<brain->mechsTypesCount; m++)
        //initialize(NrnThread nt, Memb_list * ml, type);

    //AFTER_INITIAL
    //for (int m=0; m<brain->mechsTypesCount; m++)
        //beforeAfterFunctions[AFTER_INITIAL](NrnThread nt, Memb_list * ml, type);

    //TODO nrn_deliver_events()

    //TODO nrn_deliver_events()

    //TODO nrn_spike_exchange()

    //unpin and return success
    hpx_gas_unpin(brain_addr);
    return HPX_SUCCESS;
}

hpx_action_t Brain::init = 0;
int Brain::init_handler(const int neuronsCount,
           const hpx_t neuronsAddr, const Mechanism * mechanisms,
           const size_t mechanismsCount, const int * mechDependencies)
{
    //Make sure message arrived correctly, and pin memory
    hpx_t target = hpx_thread_current_target();
    Brain *local = NULL;
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
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  init, init_handler, HPX_INT, HPX_ADDR, HPX_POINTER, HPX_SIZE_T, HPX_POINTER);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  finitialize, finitialize_handler);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  clear, clear_handler);
}

