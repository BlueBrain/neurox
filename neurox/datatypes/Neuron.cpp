#include "neurox/Neurox.h"
#include "string.h"

using namespace Neurox;

Neuron::~Neuron()
{
    delete [] synapses;
}

hpx_action_t Neuron::finitialize=0;
int Neuron::finitialize_handler()
{
    neurox_hpx_pin(Neuron);
    //set up by finitialize.c:nrn_finitialize() -> fadvance_core.c:dt2thread()
    local->cj = inputParams->secondorder ? 2.0/inputParams->dt : 1.0/inputParams->dt;

    //set up by finitialize.c:nrn_finitialize(): if (setv)
    double v = inputParams->voltage;
    hpx_call_sync(local->soma, Branch::setV, NULL, 0, v);

    // the INITIAL blocks are ordered so that mechanisms that write
    // concentrations are after ions and before mechanisms that read
    // concentrations.
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::modFunctionId::before_initialize);
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::modFunctionId::initialize);
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::modFunctionId::after_initialize);

    Neuron::setupTreeMatrixMinimal(local);

    neurox_hpx_unpin;
}

void Neuron::setupTreeMatrixMinimal(Neuron * local)
{
    hpx_call_sync(local->soma, Branch::setupMatrixInitValues, NULL, 0);

    //finitialize.c:nrn_finitialize()->set_tree_matrix_minimal->nrn_rhs
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::modFunctionId::before_breakpoint);

    //note that CAP has no current
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::modFunctionId::current);

    //finitialize.c:nrn_finitialize()->set_tree_matrix_minimal->nrn_lhs (treeset_core.c)
    // now the internal axial currents.
    //The extracellular mechanism contribution is already done.
    //	rhs += ai_j*(vi_j - vi)
    char isSoma=1;
    double dummyVoltage=-1;
    hpx_call_sync(local->soma, Branch::setupMatrixRHS, NULL, 0, isSoma, dummyVoltage);

    // calculate left hand side of
    //cm*dvm/dt = -i(vm) + is(vi) + ai_j*(vi_j - vi)
    //cx*dvx/dt - cm*dvm/dt = -gx*(vx - ex) + i(vm) + ax_j*(vx_j - vx)
    //with a matrix so that the solution is of the form [dvm+dvx,dvx] on the right
    //hand side after solving.
    //This is a common operation for fixed step, cvode, and daspk methods
    // note that CAP has no jacob
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::modFunctionId::jacob);

    //finitialize.c:nrn_finitialize()->set_tree_matrix_minimal->nrn_rhs (treeset_core.c)
    //now the cap current can be computed because any change to cm
    //by another model has taken effect. note, the first is CAP
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::modFunctionId::nrn_cap_jacob);

    //now add the axial currents
    hpx_call_sync(local->soma, Branch::setupMatrixLHS, NULL, 0, isSoma);
}

hpx_action_t Neuron::init = 0;
int Neuron::init_handler(const int gid, const hpx_t topBranch,
    double APThreshold)
{
    neurox_hpx_pin(Neuron);
    //copy information
    local->t=0; //TODO should be from InputParams
    local->topBranch=topBranch;
    local->id=gid;
    local->APThreshold=APThreshold;
    neurox_hpx_unpin;
}

hpx_action_t Neuron::addSynapses = 0;
int Neuron::addSynapses_handler(Synapse * synapses, size_t synapsesCount)
{
    neurox_hpx_pin(Neuron);
    local->synapsesCount=synapsesCount;
    local->synapses = new Synapse[synapsesCount];
    memcpy(local->synapses, synapses, sizeof(Synapse)*synapsesCount);
    neurox_hpx_unpin;
}

hpx_t Neuron::fireActionPotential(Neuron * local)
{
    //netcvode.cpp::PreSyn::send()
    Synapse *& synapses = local->synapses;
    hpx_t spikesLco = hpx_lco_and_new(local->synapsesCount);
    hpx_par_for( [&] (int i, void*)
    {
        synapses[i].deliveryTime = local->t+synapses[i].delay;
        hpx_call(local->synapses[i], Branch::queueSpike,
                 spikesLco, local->id, sizeof(int));
    },
    0, local->synapsesCount, NULL);
    return spikesLco;
}


void Neuron::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED, init, init_handler, HPX_INT, HPX_ADDR, HPX_DOUBLE, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED, addSynapses, addSynapses_handler, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE, finitialize, finitialize_handler);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE, solve, solve_handler);
}
