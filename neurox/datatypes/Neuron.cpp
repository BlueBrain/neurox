#include "neurox/Neurox.h"
#include "string.h"

using namespace NeuroX;
using namespace NeuroX::Solver;

Neuron::~Neuron()
{
    synapsesMutex = hpx_lco_sema_new(1);
}

double Neuron::getSomaVoltage()
{
    //TODO make local call
    double Vm;
    hpx_call_sync(soma, Branch::getSomaVoltage, &Vm, sizeof(double));
    return Vm;
}

void Neuron::callModFunction(Mechanism::ModFunction functionId)
{
    //TODO is it worth to make a local call instead of a call to hpx address?
    //same will all calls to soma hpx_t eg. callModFunction_handler
    assert(this==DEBUG_NEURON_DELETE);
    hpx_call_sync(soma, Branch::callModFunction, NULL, 0, &functionId, sizeof(functionId));
}

void Neuron::callNetReceiveFunction(char isInitFunction)
{
    //TODO fix
    //hpx_call_sync(soma, Branch::callNetReceiveFunction, NULL, 0,
                  //&functionId, sizeof(functionId), &t, sizeof(t), &dt, sizeof(dt);
}

hpx_action_t Neuron::finitialize=0;
int Neuron::finitialize_handler()
{
    neurox_hpx_pin(Neuron);
    assert(local==DEBUG_NEURON_DELETE);
    //set up by finitialize.c:nrn_finitialize() -> fadvance_core.c:dt2thread()
    local->cj = inputParams->secondorder ? 2.0/inputParams->dt : 1.0/inputParams->dt;

    //set up by finitialize.c:nrn_finitialize(): if (setv)
    double v = inputParams->voltage;
    hpx_call_sync(local->soma, Branch::setV, NULL, 0, &v, sizeof(v));

    // the INITIAL blocks are ordered so that mechanisms that write
    // concentrations are after ions and before mechanisms that read
    // concentrations.
    //TODO we can call this at branch level, run the three inside branch, not sequentially at neuron
    local->callModFunction(Mechanism::ModFunction::before_initialize);
    local->callModFunction(Mechanism::ModFunction::initialize);
    local->callModFunction(Mechanism::ModFunction::after_initialize);
    local->setupTreeMatrixMinimal();

    neurox_hpx_unpin;
}

void Neuron::setupTreeMatrixMinimal()
{
    hpx_call_sync(soma, Branch::setupMatrixInitValues, NULL, 0);

    //finitialize.c:nrn_finitialize()->set_tree_matrix_minimal->nrn_rhs
    callModFunction(Mechanism::ModFunction::before_breakpoint);

    //note that CAP has no current
    callModFunction(Mechanism::ModFunction::current);

    //finitialize.c:nrn_finitialize()->set_tree_matrix_minimal->nrn_lhs (treeset_core.c)
    // now the internal axial currents.
    //The extracellular mechanism contribution is already done.
    //	rhs += ai_j*(vi_j - vi)
    double dummyVoltage=-1;
    hpx_call_sync(soma, HinesSolver::setupMatrixRHS, NULL, 0,
                  &dummyVoltage, sizeof(dummyVoltage));

    // calculate left hand side of
    //cm*dvm/dt = -i(vm) + is(vi) + ai_j*(vi_j - vi)
    //cx*dvx/dt - cm*dvm/dt = -gx*(vx - ex) + i(vm) + ax_j*(vx_j - vx)
    //with a matrix so that the solution is of the form [dvm+dvx,dvx] on the right
    //hand side after solving.
    //This is a common operation for fixed step, cvode, and daspk methods
    // note that CAP has no jacob
    callModFunction(Mechanism::ModFunction::jacob);

    //finitialize.c:nrn_finitialize()->set_tree_matrix_minimal->nrn_rhs (treeset_core.c)
    //now the cap current can be computed because any change to cm
    //by another model has taken effect. note, the first is CAP
    callModFunction(Mechanism::ModFunction::capacitanceJacob);

    //now add the axial currents
    hpx_call_sync(soma, HinesSolver::setupMatrixLHS, NULL, 0);
}

hpx_action_t Neuron::init = 0;
int Neuron::init_handler(const int nargs, const void *args[], const size_t[])
{
    neurox_hpx_pin(Neuron);
    assert(nargs==3);
    DEBUG_NEURON_DELETE=local;
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

hpx_t Neuron::fireActionPotential()
{
    if (synapses.size()==0)
        return HPX_NULL;

    //netcvode.cpp::PreSyn::send()
    hpx_t spikesLco = hpx_lco_and_new(synapses.size());
    for (int s=0; s<synapses.size(); s++)
        hpx_call(synapses[s], Branch::queueSpikes, spikesLco,
                 &id, sizeof(id), &t, sizeof(t) );
    return spikesLco;
}

void Neuron::registerHpxActions()
{
    neurox_hpx_register_action(2, Neuron::init);
    neurox_hpx_register_action(1, Neuron::addSynapseTarget);
    neurox_hpx_register_action(0, Neuron::finitialize);
}
