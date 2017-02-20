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

//    bool prevFlag = synapsesTransmissionFlag;
//    synapsesTransmissionFlag = v>=threshold && synapsesTransmissionFlag==false;
//    return synapsesTransmissionFlag;
}

void Neuron::sendSpikes(spike_time_t t)
{
    //netcvode.cpp::PreSyn::send()
    if (synapsesTargets.size()>0)
    {
      //we dont have Netcon->active flag, we only add active
      //synapses to our model.
      hpx_t spikesLco = hpx_lco_and_new(synapsesTargets.size());
      for (int s=0; s<synapsesTargets.size(); s++)
      {
          //deliveryTime (t+delay) is handled on post-syn side
          hpx_call(synapsesTargets[s], Branch::addSpikeEvent, spikesLco,
                 &gid, sizeof(gid), &t, sizeof(t));
      }
      this->synapsesLCO.push_front(spikesLco);
    }
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
