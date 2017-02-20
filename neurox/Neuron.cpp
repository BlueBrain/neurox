#include "neurox/neurox.h"
#include <cstring>
#include <algorithm>

using namespace neurox;
using namespace neurox::Solver;

Neuron::Neuron(neuron_id_t neuronId, floble_t APthreshold, int thvar_index):
    gid(neuronId), threshold(APthreshold), thvar_index(thvar_index)
{
    this->synapsesMutex = hpx_lco_sema_new(1);
    this->synapsesLCO = std::deque<hpx_t> ();
    this->synapses = std::vector<hpx_t> ();
    this->synapsesTransmissionFlag = false;
}

Neuron::~Neuron() {
    hpx_lco_delete_sync(synapsesMutex);
}

size_t Neuron::getNetConsCount()
{
    hpx_lco_sema_p(synapsesMutex);
    return synapses.size();
    hpx_lco_sema_v_sync(synapsesMutex);
}

void Neuron::addSynapseTarget(hpx_t target)
{
    hpx_lco_sema_p(synapsesMutex);
    if (std::find(synapses.begin(), synapses.end(), target) == synapses.end())
    {
        synapses.push_back(target);
        synapses.shrink_to_fit();
    }
    else
    { assert(0);} //should be filtered by the branch
    hpx_lco_sema_v_sync(synapsesMutex);
}

bool Neuron::checkAPthreshold(floble_t v)
{
    //can only spike if AP threshold has been reach and I havent transmitted them already
    synapsesTransmissionFlag = v>=threshold && synapsesTransmissionFlag ==false;
    return synapsesTransmissionFlag;
}

void Neuron::sendSpikes(spike_time_t t)
{
    //netcvode.cpp::PreSyn::send()
    if (synapses.size()>0)
    {
      //we dont have Netcon->active flag, we only add active
      //synapses to our model.
      hpx_t spikesLco = hpx_lco_and_new(synapses.size());
      for (int s=0; s<synapses.size(); s++)
      {
          //deliveryTime (t+delay) is handled on post-syn side
          hpx_call(synapses[s], Branch::addSpikeEvent, spikesLco,
                 &gid, sizeof(gid), &t, sizeof(t));
      }
      this->synapsesLCO.push_front(spikesLco);
    }
    else
      this->synapsesLCO.push_front(HPX_NULL);
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
