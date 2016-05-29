#include "neurox/neurox.h"
#include "string.h"

Neuron::~Neuron()
{
    delete [] synapses;
}

hpx_action_t Neuron::init = 0;
int Neuron::init_handler(const int gid, const hpx_t topBranch,
    double APThreshold, hpx_t * synapses, size_t synapsesCount)
{
    neurox_hpx_pin(Neuron);
    //copy information
    local->t=0; //TODO should be from InputParams
    local->topBranch=topBranch;
    local->id=gid;
    local->synapsesCount=synapsesCount;
    local->APThreshold=APThreshold;
    local->synapses = new Synapse[synapsesCount];
    memcpy(local->synapses, synapses, sizeof(Synapse)*synapsesCount);
    neurox_hpx_unpin;
}

hpx_action_t Neuron::setV = 0;
int Neuron::setV_handler(const double v)
{
    neurox_hpx_pin(Neuron);
    hpx_call_sync(local->topBranch, Branch::setV, NULL, 0, v);
    neurox_hpx_unpin;
}


hpx_action_t Neuron::callMechsFunction = 0;
int Neuron::callMechsFunction_handler(const Mechanism::Function functionId)
{
    neurox_hpx_pin(Neuron);
    hpx_call_sync(local->topBranch, Branch::callMechsFunction, NULL, 0, functionId);
    neurox_hpx_unpin;
}

hpx_action_t Neuron::setupMatrixRHS = 0;
int Neuron::setupMatrixRHS_handler()
{
    neurox_hpx_pin(Neuron);
    char isSoma=1;
    double dummyVoltage=-1;
    hpx_call_sync(local->topBranch, Branch::setupMatrixRHS, NULL, 0, isSoma, dummyVoltage);
    neurox_hpx_unpin;
}

hpx_action_t Neuron::setupMatrixLHS = 0;
int Neuron::setupMatrixLHS_handler()
{
    neurox_hpx_pin(Neuron);
    char isSoma=1;
    hpx_call_sync(local->topBranch, Branch::setupMatrixLHS, NULL, 0, isSoma);
    neurox_hpx_unpin;
}

hpx_action_t Neuron::setCj = 0;
int Neuron::setCj(const double cj)
{
    neurox_hpx_pin(Neuron);
    //set up by finitialize.c:nrn_finitialize() -> fadvance_core.c:dt2thread()
    local->cj = cj;
    neurox_hpx_unpin;
}

void Neuron::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,  init, init_handler, HPX_INT, HPX_ADDR, HPX_DOUBLE, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,  setupMatrixRHS, setupMatrixRHS_handler);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,  setupMatrixLHS, setupMatrixLHS_handler);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,  setV, setV_handler, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,  setCj, setCj_handler, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,  callMechsFunction, callMechsFunction_handler, HPX_INT);
}
