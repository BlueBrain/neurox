#include "neurox/neurox.h"
#include <cstring>
#include <algorithm>
#include <utility>

using namespace neurox;
using namespace neurox::Solver;

Neuron::Neuron(neuron_id_t neuronId,
               floble_t APthreshold, int thvar_index):
    gid(neuronId), threshold(APthreshold), thvar_index(thvar_index)
{
    this->synapsesTransmissionFlag = false;
    this->synapsesMutex = hpx_lco_sema_new(1);
    this->refractoryPeriod=0;
    this->timeDependencies = inputParams->algorithm == Algorithm::BackwardEulerDebug
            ? NULL : new TimeDependencies();
}

Neuron::~Neuron()
{
    hpx_lco_delete_sync(synapsesMutex);
    delete timeDependencies;
}

Neuron::Synapse::Synapse(hpx_t addr, floble_t nextTime, floble_t delay)
    :addr(addr), nextNotificationTime(nextTime), minDelay(delay){}

size_t Neuron::getSynapseCount()
{
    size_t size = -1;
    hpx_lco_sema_p(synapsesMutex);
    size = synapses.size();
    hpx_lco_sema_v_sync(synapsesMutex);
    return size;
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
    const double teps = 1e-6;
    if (synapses.size()>0)
      for (Synapse & s : synapses)
        //deliveryTime (t+delay) is handled on post-syn side
        if (inputParams->algorithm==neurox::Algorithm::BackwardEulerDebug)
          hpx_call_sync(s.addr, Branch::addSpikeEvent, NULL, 0,
              &this->gid, sizeof(neuron_id_t), &t, sizeof(spike_time_t));
        else
        {
            spike_time_t nextNotifTime = (spike_time_t) t+s.minDelay+refractoryPeriod-teps;
            s.nextNotificationTime = (floble_t ) nextNotifTime;
            hpx_call(s.addr, Branch::addSpikeEvent, HPX_NULL,
                &this->gid, sizeof(neuron_id_t), &t, sizeof(spike_time_t),
                &nextNotifTime, sizeof(spike_time_t));
        }
}

void Neuron::sendSteppingNotification(floble_t t)
{
    //inform all dependants that need to be notified in this step (t is the 'end of step time')
    const double teps = 1e-6;
    for (Synapse & s : synapses)
        if (s.nextNotificationTime<=t) //time step just finished
        {
            assert(s.nextNotificationTime >= t-inputParams->dt);
            //next time allowed by post-syn neuron is also the time I have to notify him
            spike_time_t nextNotifTime = (spike_time_t)  t+s.minDelay-teps;
            s.nextNotificationTime = (floble_t) nextNotifTime;
#ifdef PRINT_EVENT
            printf("Notification of step: gid %d at time %.3f: informs %llu of next notif time %.3f\n",
                   this->gid, t, s.addr, s.addr, s.nextNotificationTime);
#endif

            //tell post-syn neuron how long he can proceed to until he has to wait
            hpx_call(s.addr, Branch::updateTimeDependencyValue, HPX_NULL,
                 &this->gid, sizeof(neuron_id_t), &nextNotifTime, sizeof(spike_time_t));
        }
}

Neuron::TimeDependencies::TimeDependencies()
{
    this->dependenciesMutex = hpx_lco_sema_new(1);
    this->dependenciesMinTimeCached = -1; //to be set by Neuron::updateTimeDependenciesMinTimeCached;
    this->neuronWaitingLco = hpx_lco_and_new(1); //wait for one update at a time
    this->neuronWaitingFlag = false;
}

Neuron::TimeDependencies::~TimeDependencies() {
    hpx_lco_delete_sync(dependenciesMutex);
    hpx_lco_delete_sync(neuronWaitingLco);
}

void Neuron::TimeDependencies::updateDependenciesMinTimeCached()
{
    if (dependenciesVector.size()==0) return;
    dependenciesMinTimeCached = (*std::min_element(
        dependenciesVector.begin(), dependenciesVector.end(),
        [] (pair<neuron_id_t, floble_t> const& lhs, pair<neuron_id_t, floble_t> const& rhs)
            {return lhs.second < rhs.second;} )).second;
}

void Neuron::TimeDependencies::updateTimeDependency(
        neuron_id_t srcGid, floble_t maxTimeAllowed, bool initializationPhase)
{
    hpx_lco_sema_p(dependenciesMutex);
    assert(maxTimeAllowed > dependenciesMinTimeCached);

    if (initializationPhase)
    {
        assert(dependenciesMap.find(srcGid) == dependenciesMap.end());
        dependenciesVector.push_back(std::make_pair(srcGid, maxTimeAllowed)); //time-value
        dependenciesMap[srcGid] = &(dependenciesVector.back().second); //pointer to time-value
    }
    else
    {
        assert(dependenciesMap.find(srcGid)!= dependenciesMap.end());
        if (*dependenciesMap.at(srcGid) < maxTimeAllowed) //order of msgs is not guaranteed so take only last update (highest time value)
            *dependenciesMap.at(srcGid) = maxTimeAllowed; //set time-value in vector
    }

    if (neuronWaitingFlag) //if neuron is waiting, tell him there was an update
        hpx_lco_set(neuronWaitingLco, 0, NULL, HPX_NULL, HPX_NULL);

    hpx_lco_sema_v_sync(dependenciesMutex);
}

void Neuron::TimeDependencies::waitForTimeDependencyNeurons(floble_t waitUntilTime)
{
    //if I have no dependencies... I'm free to go!
    if (dependenciesVector.size()==0) return;

    hpx_lco_sema_p(dependenciesMutex);
    assert(dependenciesMinTimeCached>0);

    if (waitUntilTime > dependenciesMinTimeCached) //if last min-value cached does not cover my full step
        updateDependenciesMinTimeCached();         //recompute min value from existing table

    while (waitUntilTime > dependenciesMinTimeCached) //if I still cant proceed
    {
        if (!neuronWaitingFlag) neuronWaitingFlag = true;
        hpx_lco_sema_v_sync(dependenciesMutex); //unlock to allow thread to run 'updateTimeDependencyTime'
        hpx_lco_wait_reset(this->neuronWaitingLco); //sleep and wait for thread to wake me up
        hpx_lco_sema_p(dependenciesMutex);
        updateDependenciesMinTimeCached(); //recompute min value from existing table
    }
    neuronWaitingFlag = false;
    hpx_lco_sema_v_sync(dependenciesMutex);
}
