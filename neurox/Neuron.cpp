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

Neuron::Synapse::Synapse(hpx_t addr, floble_t nextTime, floble_t delay, int destinationGid)
    :addr(addr), nextNotificationTime(nextTime), minDelay(delay)
{
#ifndef NDEBUG
    this->destinationGid=destinationGid;
#endif
}

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
#if !defined(NDEBUG) && defined(PRINT_EVENT)
            printf("Neuron::sendSpikes: gid %d at time %.3f, informs gid %d (%llu) of next notif time %.3f\n",
                   this->gid, t, s.destinationGid, s.addr, s.addr, s.nextNotificationTime);
#endif
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
#if !defined(NDEBUG) && defined(PRINT_EVENT)
            printf("Neuron::sendSteppingNotification: gid %d at time %.3f, informs gid %d (%llu) of next notif time %.3f\n",
                   this->gid, t, s.destinationGid, s.addr, s.addr, s.nextNotificationTime);
#endif
            //tell post-syn neuron how long he can proceed to until he has to wait
            hpx_call(s.addr, Branch::updateTimeDependencyValue, HPX_NULL,
                 &this->gid, sizeof(neuron_id_t), &nextNotifTime, sizeof(spike_time_t));
        }
}

Neuron::TimeDependencies::TimeDependencies()
{
    this->dependenciesMinTimeCached = -1; //to be set by Neuron::updateTimeDependenciesMinTimeCached;
    this->neuronWaitingFlag = false;
    libhpx_cond_init(&this->dependenciesWaitCondition);
    libhpx_mutex_init(&this->dependenciesLock);
}

Neuron::TimeDependencies::~TimeDependencies() {
    libhpx_cond_destroy(&this->dependenciesWaitCondition);
    libhpx_mutex_destroy(&this->dependenciesLock);
}

size_t Neuron::TimeDependencies::getDependenciesCount()
{
    size_t size = -1;
    libhpx_mutex_lock(&this->dependenciesLock);
    size = dependenciesMap.size();
    libhpx_mutex_unlock(&this->dependenciesLock);
    return size;
}

void Neuron::TimeDependencies::updateDependenciesMinTimeCached()
{
    assert(dependenciesMap.size()>0);
    floble_t minValue = inputParams->tstop+99999;
    for (auto & key_value : dependenciesMap)
        if (minValue < key_value.second)
            minValue = key_value.second;
    dependenciesMinTimeCached = minValue;
    /* //TODO most efficient way to calculate this!!
    dependenciesMinTimeCached = *std::min_element(
        dependenciesMap.begin(), dependenciesMap.end(),
        [] (floble_t const& lhs, floble_t const& rhs)
            {return lhs < rhs;} );*/
}

void Neuron::TimeDependencies::updateTimeDependency(
        neuron_id_t srcGid, floble_t dependencyNotificationTime, bool initializationPhase)
{
    libhpx_mutex_lock(&this->dependenciesLock);
    if (initializationPhase)
    {
        assert(dependenciesMap.find(srcGid) == dependenciesMap.end());
        dependenciesMap[srcGid] = dependencyNotificationTime;
    }
    else
    {
        assert(dependenciesMap.find(srcGid) != dependenciesMap.end());
        if (dependenciesMap.at(srcGid) < dependencyNotificationTime) //order of msgs is not guaranteed so take only last update (highest time value)
            dependenciesMap.at(srcGid) = dependencyNotificationTime; //set time-value in vector

        if (neuronWaitingFlag) //if neuron is waiting, tell him there was an update
            libhpx_cond_broadcast(&this->dependenciesWaitCondition);
    }
    libhpx_mutex_unlock(&this->dependenciesLock);
}

void Neuron::TimeDependencies::waitForTimeDependencyNeurons(floble_t waitUntilTime, int gid)
{
    //if I have no dependencies... I'm free to go!
    if (dependenciesMap.size()==0) return;

    libhpx_mutex_lock(&this->dependenciesLock);
    assert(dependenciesMinTimeCached>0);

    if (waitUntilTime > dependenciesMinTimeCached) //if last min-value cached does not cover my full step
        updateDependenciesMinTimeCached();         //recompute min value from existing table

    while (waitUntilTime > dependenciesMinTimeCached) //if I still cant proceed
    {
        if (!neuronWaitingFlag) neuronWaitingFlag = true;

#if !defined(NDEBUG) && defined(PRINT_EVENT)
        printf("Neuron::waitForTimeDependencyNeurons: gid %d WAITING to update t>%.8f (dependenciesMinTimeCached=%.8f)\n",
               gid, waitUntilTime, dependenciesMinTimeCached);
#endif

        libhpx_cond_wait(&this->dependenciesWaitCondition, &this->dependenciesLock);
        updateDependenciesMinTimeCached(); //recompute min value from existing table

#if !defined(NDEBUG) && defined(PRINT_EVENT)
        printf("Neuron::waitForTimeDependencyNeurons: gid %d woke up (dependenciesMinTimeCached=%.8f)\n",
               gid, dependenciesMinTimeCached);
#endif
    }
    neuronWaitingFlag = false;
    libhpx_mutex_unlock(&this->dependenciesLock);
}
