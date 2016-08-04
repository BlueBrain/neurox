#include "neurox/Neurox.h"
#include "string.h"

using namespace Neurox;

Neuron::~Neuron()
{
    synapsesMutex = hpx_lco_sema_new(1);
}

double Neuron::getSomaVoltage()
{
    double Vm;
    hpx_call_sync(soma, Branch::getSomaVoltage, &Vm, sizeof(double));
    return Vm;
}

void Neuron::callModFunction(Mechanism::ModFunction functionId)
{
    //TODO fix
     //hpx_call_sync(soma, Branch::callModFunction, NULL, 0, &functionId, sizeof(functionId));
}

void Neuron::callNetReceiveFunction(char isInitFunction)
{
    //TODO fix
     //hpx_call_sync(soma, Branch::callModFunction, NULL, 0, &functionId, sizeof(functionId), t, dt);
}

hpx_action_t Neuron::finitialize=0;
int Neuron::finitialize_handler()
{
    neurox_hpx_pin(Neuron);
    //set up by finitialize.c:nrn_finitialize() -> fadvance_core.c:dt2thread()
    local->cj = inputParams->secondorder ? 2.0/inputParams->dt : 1.0/inputParams->dt;

    //set up by finitialize.c:nrn_finitialize(): if (setv)
    double v = inputParams->voltage;
    hpx_call_sync(local->soma, Branch::setV, NULL, 0, &v, sizeof(v));

    // the INITIAL blocks are ordered so that mechanisms that write
    // concentrations are after ions and before mechanisms that read
    // concentrations.
    local->callModFunction(Mechanism::ModFunction::before_initialize);
    local->callModFunction(Mechanism::ModFunction::initialize);
    local->callModFunction(Mechanism::ModFunction::after_initialize);
    Neuron::setupTreeMatrixMinimal(local);

    neurox_hpx_unpin;
}

void Neuron::setupTreeMatrixMinimal(Neuron * local)
{
    hpx_call_sync(local->soma, Branch::setupMatrixInitValues, NULL, 0);

    //finitialize.c:nrn_finitialize()->set_tree_matrix_minimal->nrn_rhs
    local->callModFunction(Mechanism::ModFunction::before_breakpoint);

    //note that CAP has no current
    local->callModFunction(Mechanism::ModFunction::current);

    //finitialize.c:nrn_finitialize()->set_tree_matrix_minimal->nrn_lhs (treeset_core.c)
    // now the internal axial currents.
    //The extracellular mechanism contribution is already done.
    //	rhs += ai_j*(vi_j - vi)
    char isSoma=1;
    double dummyVoltage=-1;
    hpx_call_sync(local->soma, Branch::setupMatrixRHS, NULL, 0,
                  &isSoma, sizeof(isSoma), &dummyVoltage, sizeof(dummyVoltage));

    // calculate left hand side of
    //cm*dvm/dt = -i(vm) + is(vi) + ai_j*(vi_j - vi)
    //cx*dvx/dt - cm*dvm/dt = -gx*(vx - ex) + i(vm) + ax_j*(vx_j - vx)
    //with a matrix so that the solution is of the form [dvm+dvx,dvx] on the right
    //hand side after solving.
    //This is a common operation for fixed step, cvode, and daspk methods
    // note that CAP has no jacob
    local->callModFunction(Mechanism::ModFunction::jacob);

    //finitialize.c:nrn_finitialize()->set_tree_matrix_minimal->nrn_rhs (treeset_core.c)
    //now the cap current can be computed because any change to cm
    //by another model has taken effect. note, the first is CAP
    local->callModFunction(Mechanism::ModFunction::capJacob);

    //now add the axial currents
    hpx_call_sync(local->soma, Branch::setupMatrixLHS, NULL, 0, &isSoma, sizeof(isSoma));
}

hpx_action_t Neuron::init = 0;
int Neuron::init_handler(const int nargs, const void *args[], const size_t[])
{
    neurox_hpx_pin(Neuron);
    assert(nargs==3);
    const int neuronId = *(const int*) args[0];
    const hpx_t topBranch = *(const hpx_t*) args[1];
    const double APthreshold = *(const double*) args[2];

    local->t = inputParams->t;
    local->dt = inputParams->dt;
    local->soma = topBranch;
    local->id = neuronId;
    local->APthreshold = APthreshold;
    neurox_hpx_unpin;
}

hpx_action_t Neuron::addSynapseTarget = 0;
int Neuron::addSynapseTarget_handler(const hpx_t * synapseTarget, const size_t)
{
    neurox_hpx_pin(Neuron);
    hpx_lco_sema_p(local->synapsesMutex);
    local->synapses.push_back(*synapseTarget);
    local->synapses.shrink_to_fit();
    hpx_lco_sema_v_sync(local->synapsesMutex);
    neurox_hpx_unpin;
}

hpx_t Neuron::fireActionPotential(Neuron * local)
{
    if (local->synapses.size()==0)
        return HPX_NULL;

    //netcvode.cpp::PreSyn::send()
    hpx_t spikesLco = hpx_lco_and_new(local->synapses.size());
    for (int s=0; s<local->synapses.size(); s++)
        hpx_call(local->synapses[s], Branch::queueSpikes, spikesLco,
                 &local->id, sizeof(local->id), &local->t, sizeof(local->t) );
    return spikesLco;
}

void Neuron::registerHpxActions()
{
    neurox_hpx_register_action(2, Neuron::init);
    neurox_hpx_register_action(1, Neuron::addSynapseTarget);
    neurox_hpx_register_action(0, Neuron::finitialize);
}
