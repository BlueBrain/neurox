#include "neurox/neurox.h"
#include <cstring>
#include <algorithm>

using namespace neurox;
using namespace neurox::Solver;

Neuron::Neuron(neuron_id_t neuronId,
               floble_t APthreshold, int thvar_index):
    gid(neuronId), threshold(APthreshold), thvar_index(thvar_index)
{
    this->synapsesMutex = hpx_lco_sema_new(1);
    this->synapsesTargets = std::vector<hpx_t> ();
    this->synapsesTransmissionFlag = false;
    this->slidingTimeWindow = new SlidingTimeWindow(
        MIN_DELAY_IN_INTERVALS*inputParams->dt);
}

Neuron::~Neuron() {
    hpx_lco_delete_sync(synapsesMutex);
}

size_t Neuron::getSynapsesCount()
{
    hpx_lco_sema_p(synapsesMutex);
    return synapsesTargets.size();
    hpx_lco_sema_v_sync(synapsesMutex);
}

void Neuron::addSynapseTarget(hpx_t target)
{
    hpx_lco_sema_p(synapsesMutex);
    if (std::find(synapsesTargets.begin(), synapsesTargets.end(), target) == synapsesTargets.end())
    {
        synapsesTargets.push_back(target);
        synapsesTargets.shrink_to_fit();
    }
    else
    { assert(0);} //should be filtered by the branch
    hpx_lco_sema_v_sync(synapsesMutex);
}

//netcvode.cpp::static bool pscheck(...)
bool Neuron::checkAPthresholdAndTransmissionFlag(floble_t v)
{
    //can only spike if AP threshold has been reach and spikes havent already been transmitted
    if (v > threshold) {
        if (synapsesTransmissionFlag == false) {
            synapsesTransmissionFlag = true;
            return true;
        }
    } else {
        synapsesTransmissionFlag = false;
    }
    return false;
}

void Neuron::sendSpikes(spike_time_t t)
{
    //netcvode.cpp::PreSyn::send()
    if (synapsesTargets.size()>0)
      for (hpx_t & destinationAddr : synapsesTargets)
        //deliveryTime (t+delay) is handled on post-syn side
        hpx_call(destinationAddr, Branch::addSpikeEvent, HPX_NULL,
                 &gid, sizeof(gid), &t, sizeof(t));

        //hpx_call_sync(destinationAddr, Branch::addSpikeEvent,
        //       NULL, 0, &gid, sizeof(gid), &t, sizeof(t));
}

Neuron::SlidingTimeWindow::SlidingTimeWindow (floble_t windowSize)
    :windowSize(windowSize)
{
    //initial time of all dependencis
    dependencies_t = inputParams->tstart;

    //depedent_lco is the hpx_addr of someone else's dependancy_lco
    //(others will share with me)
    dependant_lco = HPX_NULL;

    //dependancy_lco is stored locally so dependency neurons write in it;
    //(to be shared with others)
    dependencies_lco = hpx_lco_and_new(1);
}

Neuron::SlidingTimeWindow::~SlidingTimeWindow()
{} //TODO

void Neuron::SlidingTimeWindow::informTimeDependants(floble_t t)
{
    //inform whoever is dependant on my that i started a new step
    static const double teps = 1e-10;
    floble_t time = t+teps;
    hpx_lco_set(dependant_lco, sizeof(time), &time, HPX_NULL, HPX_NULL);
}

void Neuron::SlidingTimeWindow::waitForTimeDependencies(floble_t t)
{
    //wait for my dependency to be at least a 'min delay' behind me
    while (dependencies_t < t-windowSize);
        hpx_lco_get_reset(dependencies_lco , sizeof(floble_t), &dependencies_t);
}

hpx_action_t Neuron::SlidingTimeWindow::initDependencies = 0;
int Neuron::SlidingTimeWindow::initDependencies_handler()
{
    neurox_hpx_pin(Branch);
    assert(local->soma);

    //share HPX resources
    int nrnThreadId = local->nt->id;
    hpx_t addr = getNeuronAddr(nrnThreadId==0 ? neurox::neuronsCount-1 : nrnThreadId-1);
    hpx_call_sync(addr, Neuron::SlidingTimeWindow::setDependant, HPX_NULL, 0,
                  &local->soma->slidingTimeWindow->dependencies_lco, sizeof(hpx_t));

    neurox_hpx_unpin;
}

hpx_action_t Neuron::SlidingTimeWindow::setDependant = 0;
int Neuron::SlidingTimeWindow::setDependant_handler(const hpx_t * dependant_ptr, const size_t)
{
    neurox_hpx_pin(Branch);
    assert(local->soma);
    local->soma->slidingTimeWindow->dependant_lco = *dependant_ptr;
    neurox_hpx_unpin;
}

void Neuron::registerHpxActions()
{
    neurox_hpx_register_action(0, Neuron::SlidingTimeWindow::initDependencies);
    neurox_hpx_register_action(1, Neuron::SlidingTimeWindow::setDependant);
}
