#include "neurox/Neurox.h"
#include "string.h"

using namespace Neurox;

Neuron::~Neuron()
{
    synapsesMutex = hpx_lco_sema_new(1);
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
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::Functions::before_initialize, local->t, local->dt);
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::Functions::initialize, local->t, local->dt);
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::Functions::after_initialize, local->t, local->dt);

    Neuron::setupTreeMatrixMinimal(local);

    neurox_hpx_unpin;
}

void Neuron::setupTreeMatrixMinimal(Neuron * local)
{
    hpx_call_sync(local->soma, Branch::setupMatrixInitValues, NULL, 0);

    //finitialize.c:nrn_finitialize()->set_tree_matrix_minimal->nrn_rhs
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::Functions::before_breakpoint, local->t, local->dt);

    //note that CAP has no current
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::Functions::current, local->t, local->dt);

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
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::Functions::jacob, local->t, local->dt);

    //finitialize.c:nrn_finitialize()->set_tree_matrix_minimal->nrn_rhs (treeset_core.c)
    //now the cap current can be computed because any change to cm
    //by another model has taken effect. note, the first is CAP
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::Functions::capJacob, local->t, local->dt);

    //now add the axial currents
    hpx_call_sync(local->soma, Branch::setupMatrixLHS, NULL, 0, isSoma);
}

hpx_action_t Neuron::init = 0;
int Neuron::init_handler(const int gid, const hpx_t topBranch, double APthreshold)
{
    neurox_hpx_pin(Neuron);
    //copy information
    local->t=inputParams->t;
    local->dt=inputParams->dt;
    local->soma=topBranch;
    local->id=gid;
    local->APthreshold=  APthreshold;
    neurox_hpx_unpin;
}

hpx_action_t Neuron::addSynapse = 0;
int Neuron::addSynapse_handler(const hpx_t synapseTarget)
{
    neurox_hpx_pin(Neuron);
    hpx_lco_sema_p(local->synapsesMutex);
    local->synapses.push_back(synapseTarget);
    local->synapses.shrink_to_fit();
    hpx_lco_sema_v_sync(local->synapsesMutex);
    neurox_hpx_unpin;
}

hpx_t Neuron::fireActionPotential(Neuron * local)
{
    if (local->synapses.size()==0) return HPX_NULL;

    //netcvode.cpp::PreSyn::send()
    hpx_t spikesLco = hpx_lco_and_new(local->synapses.size());
    for (int s=0; s<local->synapses.size(); s++)
        hpx_call(local->synapses[s], Branch::queueSpikes,
                 spikesLco, local->id, local->t );
    return spikesLco;
}


void Neuron::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED, init, init_handler, HPX_INT, HPX_ADDR, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED, addSynapse, addSynapse_handler, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE, finitialize, finitialize_handler);
}
