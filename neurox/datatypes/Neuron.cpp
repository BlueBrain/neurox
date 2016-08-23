#include "neurox/Neurox.h"
#include <cstring>
#include <algorithm>

using namespace NeuroX;
using namespace NeuroX::Solver;

Neuron::~Neuron() {
    hpx_lco_delete_sync(synapsesMutex);
}

double Neuron::getSomaVoltage()
{
    //TODO make local call (?)
    double Vm;
    hpx_call_sync(soma, Branch::getSomaVoltage, &Vm, sizeof(double));
    return Vm;
}

void Neuron::callModFunction(Mechanism::ModFunction functionId)
{
    //TODO is it worth to make a local call instead of a call to hpx address?
    //same will all calls to soma hpx_t eg. callModFunction_handler
    if (functionId<BEFORE_AFTER_SIZE) return; //not used
    hpx_call_sync(soma, Branch::callModFunction, NULL, 0, &functionId, sizeof(functionId));
}

void Neuron::setupTreeMatrixMinimal()
{
    hpx_call_sync(soma, HinesSolver::setupMatrixRHS, NULL, 0);
    hpx_call_sync(soma, HinesSolver::gaussianFwdTriangulation, NULL, 0, NULL,0);
    hpx_call_sync(soma, HinesSolver::setupMatrixLHS, NULL, 0);
    hpx_call_sync(soma, HinesSolver::gaussianFwdTriangulation, NULL, 0, NULL,0);
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
    local->synapsesMutex = hpx_lco_sema_new(1);
    neurox_hpx_unpin;
}

hpx_action_t Neuron::addSynapseTarget = 0;
int Neuron::addSynapseTarget_handler(const hpx_t * synapseTarget, const size_t)
{
    neurox_hpx_pin(Neuron);
    hpx_lco_sema_p(local->synapsesMutex);
    if (std::find(local->synapses.begin(), local->synapses.end(), *synapseTarget) == local->synapses.end())
    {
        local->synapses.push_back(*synapseTarget);
        local->synapses.shrink_to_fit();
    }
    else
    { assert(0);} //should be filtered by the branch
    hpx_lco_sema_v_sync(local->synapsesMutex);
    neurox_hpx_unpin;
}

void Neuron::fireActionPotential()
{
    //netcvode.cpp::PreSyn::send()
    hpx_t spikesLco = hpx_lco_and_new(synapses.size());
    for (int s=0; s<synapses.size(); s++)
    {
        double tt = this->t +  1e-10;
        hpx_call(synapses[s], Branch::queueSpikes, spikesLco,
                 &id, sizeof(id), &tt, sizeof(tt) );
    }
    this->synapsesLCO.push_front(spikesLco);
}

void Neuron::waitForSynapsesDelivery(int commStepSize)
{
    assert(this->synapsesLCO.size()<=commStepSize);
    if (this->synapsesLCO.size()==commStepSize)
    {
        hpx_lco_wait(this->synapsesLCO.back());
        this->synapsesLCO.pop_back();
    }
}

hpx_action_t Neuron::broadcastNetCons = 0;
int Neuron::broadcastNetCons_handler(const int nargs, const void *args[], const size_t sizes[])
{
    neurox_hpx_pin(Neuron);
    assert(nargs==2);
    hpx_call_sync(local->soma, Branch::broadcastNetCons, NULL, 0,
                  (void*) args[0], sizes[0], (void*) args[1], sizes[1]);
    neurox_hpx_unpin;
}

void Neuron::registerHpxActions()
{
    neurox_hpx_register_action(2, Neuron::init);
    neurox_hpx_register_action(2, Neuron::broadcastNetCons);
    neurox_hpx_register_action(1, Neuron::addSynapseTarget);
}
