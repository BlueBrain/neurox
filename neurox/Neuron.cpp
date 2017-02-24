#include "neurox/neurox.h"
#include <cstring>
#include <algorithm>

using namespace neurox;
using namespace neurox::Solver;

Neuron::Neuron(neuron_id_t neuronId, floble_t APthreshold, int thvar_index):
    gid(neuronId), threshold(APthreshold), thvar_index(thvar_index)
{
    this->synapsesMutex = hpx_lco_sema_new(1);
    this->slidingTimeWindow = new SlidingTimeWindow(MIN_DELAY_IN_INTERVALS);
    this->synapsesTargets = std::vector<hpx_t> ();
    this->synapsesTransmissionFlag = false;
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
    {
      hpx_t spikesLco = hpx_lco_and_new(synapsesTargets.size());
      for (hpx_t & destinationAddr : synapsesTargets)
      {
          //deliveryTime (t+delay) is handled on post-syn side
          //hpx_call(destinationAddr, Branch::addSpikeEvent, spikesLco,
          //       &gid, sizeof(gid), &t, sizeof(t));

          hpx_call_sync(destinationAddr, Branch::addSpikeEvent,
                 NULL, 0, &gid, sizeof(gid), &t, sizeof(t));
      }
      this->slidingTimeWindow->lcos.push_front(spikesLco);
    }
    //TODO this else is never hit i think!
    else
      this->slidingTimeWindow->lcos.push_front(HPX_NULL);
}


Neuron::SlidingTimeWindow::SlidingTimeWindow (size_t windowSize)
{
    for (int i=0; i<windowSize-1; i++)
        this->lcos.push_front(HPX_NULL);
}


void Neuron::SlidingTimeWindow::waitForSlidingTimeWindow()
{
    assert(lcos.size() == MIN_DELAY_IN_INTERVALS);
    hpx_t lastWindow = lcos.back();
    lcos.pop_back();

    if (lastWindow != HPX_NULL) //ie if it spiked
    {
      hpx_lco_wait(lastWindow);
      hpx_lco_delete_sync(lastWindow);
    }
}
