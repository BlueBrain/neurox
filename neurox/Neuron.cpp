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
    this->timeDependenciesMutex = hpx_lco_sema_new(1);
    this->synapsesTransmissionFlag = false;
}

Neuron::~Neuron() {
    hpx_lco_delete_sync(synapsesMutex);
    hpx_lco_delete_sync(timeDependenciesMutex);
}

Neuron::Synapse::Synapse(hpx_t addr, floble_t nextTime, floble_t delay)
    :addr(addr), nextNotificationTime(nextTime), minDelay(delay){}

size_t Neuron::getSynapsesCount()
{
    hpx_lco_sema_p(synapsesMutex);
    return synapses.size();
    hpx_lco_sema_v_sync(synapsesMutex);
}

void Neuron::addSynapse(Synapse syn)
{
    hpx_lco_sema_p(synapsesMutex);
    synapses.push_back(syn);
    synapses.shrink_to_fit();
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
    if (synapses.size()>0)
      for (Synapse & s : synapses)
        //deliveryTime (t+delay) is handled on post-syn side
        if (inputParams->algorithm==neurox::Algorithm::BackwardEulerDebug)
          hpx_call_sync(s.addr, Branch::addSpikeEvent, NULL, 0,
              &this->gid, sizeof(neuron_id_t), &t, sizeof(spike_time_t));
        else
          hpx_call(s.addr, Branch::addSpikeEvent, HPX_NULL,
              &this->gid, sizeof(neuron_id_t), &t, sizeof(spike_time_t));
}

void Neuron::informTimeDependantNeurons(floble_t t)
{
    //inform all dependants that need to be notified in this step
    static const double teps = 1e-10;
    for (Synapse & s : synapses)
        if (s.nextNotificationTime<t) //time step just finished
        {
            //next time allowed by post-syn neuron is also the time I have to notify him
            s.nextNotificationTime = t+s.minDelay-teps;

            //tell post-syn neuron how long he can proceed to until he has to wait
            hpx_call(s.addr, Branch::updateTimeDependencyTime, HPX_NULL,
                 &this->gid, sizeof(neuron_id_t), &s.nextNotificationTime, sizeof(floble_t));
        }
}

void Neuron::updateTimeDependenciesMinTimeChached()
{
    floble_t minTime = inputParams->tstop;
    hpx_lco_sema_p(timeDependenciesMutex);
    for (auto & it : timeDependencies)
        minTime = std::min(minTime,it.second);
    hpx_lco_sema_v_sync(timeDependenciesMutex);
    timeDependenciesMinTimeCached = minTime;
}

void Neuron::updateTimeDependencyTime(neuron_id_t srcGid, floble_t maxTimeAllowed)
{
    assert(maxTimeAllowed > timeDependenciesMinTimeCached);
    hpx_lco_sema_p(timeDependenciesMutex);
    assert(timeDependencies.find(srcGid)!= timeDependencies.end());
    timeDependencies.at(srcGid) = maxTimeAllowed;
    hpx_lco_sema_v_sync(timeDependenciesMutex);
}

void Neuron::waitForTimeDependencyNeurons(floble_t waitUntilTime)
{
    //if my knowledge of max-time allowed is too small for me to step,
    //iterate the latest information from dependencies to recompute it;
    while(waitUntilTime > timeDependenciesMinTimeCached)
        updateTimeDependenciesMinTimeChached();
}
