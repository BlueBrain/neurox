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

hpx_action_t Neuron::addOutgoingSynapses = 0;
int Neuron::addOutgoingSynapses_handler(hpx_t * synapsesTargets, size_t synapsesCount)
{
    neurox_hpx_pin(Neuron);
    local->synapsesCount=synapsesCount;
    local->synapsesTargets = new hpx_t[synapsesCount];
    memcpy(local->synapses, synapses, sizeof(hpx_t)*synapsesCount);
    neurox_hpx_unpin;
}

/* TODO: for the synapses transmission:
 * Neuron only sends GID and delivery time as 2 bytes (so that data fits in the header).
 * Michael Hines says that sending the struct takes too much bandwith, therefore
 * they store the struct on the POST-SYNAPTIC side, and on the pre they only send the identifier
 * (ie pre-syn gid+delivery time)
 */
hpx_t Neuron::fireActionPotential(Neuron * local)
{
    //netcvode.cpp::PreSyn::send()
    Synapse *& synapses = local->synapses;
    hpx_t spikesLco = hpx_lco_and_new(local->synapsesCount);
    for (int s=0; s<local->synapsesCount; s++)
    {
        double deliveryTime = local->t+synapses[s].delay;
        hpx_call(local->synapsesTargets[s], Branch::queueSpike,
                 spikesLco, local->id, deliveryTime );
    }
    return spikesLco;
}


void Neuron::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED, init, init_handler, HPX_INT, HPX_ADDR, HPX_DOUBLE, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED, addOutgoingSynapses, addOutgoingSynapses_handler, HPX_ADDR, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED, addIncomingSynapse, addIncomingSynapse_handle, HPX_POINTER, HPX_DOUBLE, HPX_DOUBLE, HPX_INT, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE, finitialize, finitialize_handler);
}
