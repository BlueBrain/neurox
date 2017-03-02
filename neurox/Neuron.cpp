#include "neurox/neurox.h"
#include <cstring>
#include <algorithm>
#include <utility>

using namespace neurox;
using namespace neurox::Solver;

const double teps = 1e-6;

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

Neuron::Synapse::Synapse(hpx_t addr, floble_t minDelay, int destinationGid)
    //TODO :addr(addr), minDelay(delay)
    :addr(addr), minDelay(0.05)
{
    this->nextNotificationTime=this->minDelay-teps;
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
    if (synapses.size()>0)
      for (Synapse & s : synapses)
        if (inputParams->algorithm==neurox::Algorithm::BackwardEulerDebug)
        {
            //deliveryTime (t+delay) is handled on post-syn side (diff value for every NetCon)
            hpx_call_sync(s.addr, Branch::addSpikeEvent, NULL, 0,
                &this->gid, sizeof(neuron_id_t), &t, sizeof(spike_time_t));
        }
        else
        {
            //s.minDelay is the min delay accross all Netcons on the post-syn neuron
            spike_time_t nextNotifTime = (spike_time_t) t+s.minDelay+refractoryPeriod-teps;
            //TODO s.nextNotificationTime = (floble_t ) nextNotifTime;
            hpx_call(s.addr, Branch::addSpikeEvent, HPX_NULL,
                &this->gid, sizeof(neuron_id_t), &t, sizeof(spike_time_t),
                &nextNotifTime, sizeof(spike_time_t));
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
            printf("Neuron::sendSpikes: gid %d at time %.3f, informs gid %d of next notif time %.3f+%.3f+%.3f=%.3f\n",
                   this->gid, t, s.destinationGid, t, s.minDelay, refractoryPeriod, s.nextNotificationTime);
#endif
        }
}

void Neuron::sendSteppingNotification(floble_t timeEndStep)
{
    //inform all dependants that need to be notified in this step (t is the 'end of step time')
    for (Synapse & s : synapses)
        if (s.nextNotificationTime<=timeEndStep) //time step just finished
        {
            assert(s.nextNotificationTime >= timeEndStep-inputParams->dt); //must have been covered by previous steps
            //next time allowed by post-syn neuron is also the time I have to notify him
            spike_time_t nextNotifTime = (spike_time_t)  timeEndStep+s.minDelay-teps;
            s.nextNotificationTime = (floble_t) nextNotifTime;
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
            printf("Neuron::sendSteppingNotification: gid %d at time %.3f, informs gid %d of next notif time %.3f+%.3f=%.8f\n",
                   this->gid, timeEndStep, s.destinationGid, timeEndStep, s.minDelay, s.nextNotificationTime);
#endif
            //tell post-syn neuron how long he can proceed to until he has to wait
            hpx_call(s.addr, Branch::updateTimeDependencyValue, HPX_NULL,
                 &this->gid, sizeof(neuron_id_t), &nextNotifTime, sizeof(spike_time_t));
        }
}

Neuron::TimeDependencies::TimeDependencies()
{
    libhpx_cond_init(&this->dependenciesWaitCondition);
    libhpx_mutex_init(&this->dependenciesLock);
    this->dependenciesTimeNeuronWaitsFor = 0; //0 means not waiting
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


floble_t Neuron::TimeDependencies::getDependenciesMinTime()
{
    assert(dependenciesMap.size()>0);
    return std::min_element(
        dependenciesMap.begin(), dependenciesMap.end(),
        [] (pair<neuron_id_t, floble_t> const& lhs, pair<neuron_id_t, floble_t> const& rhs)
            {return lhs.second < rhs.second;} )->second;
}

void Neuron::TimeDependencies::updateTimeDependency(
        neuron_id_t srcGid, floble_t dependencyNotificationTime, bool initializationPhase)
{

#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
    printf("Neuron::updateTimeDependency: srcGid %d waits for lock (dependencyNotificationTime=%.8f)\n",
           srcGid, dependencyNotificationTime);
#endif
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
            dependenciesMap.at(srcGid) = dependencyNotificationTime;

        if (dependenciesTimeNeuronWaitsFor > 0) //if neuron is waiting for a dependencies time update
          if (getDependenciesMinTime() >= dependenciesTimeNeuronWaitsFor) //and new min time allows neuron to proceed
          {
            dependenciesTimeNeuronWaitsFor = 0; //mark neuron as not asleep anymore
            libhpx_cond_broadcast(&this->dependenciesWaitCondition); //wake up neuron
          }
    }
    libhpx_mutex_unlock(&this->dependenciesLock);
}

void Neuron::TimeDependencies::waitForTimeDependencyNeurons(floble_t endOfStepTime, int gid)
{
    //if I have no dependencies... I'm free to go!
    if (dependenciesMap.size()==0) return;

#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
    printf("Neuron::waitForTimeDependencyNeurons: gid %d STARTS (t+dt=%.8f)\n", gid, endOfStepTime);
#endif
    libhpx_mutex_lock(&this->dependenciesLock);

    if (getDependenciesMinTime() < endOfStepTime) //if I cant proceed
    {
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
        printf("Neuron::waitForTimeDependencyNeurons: gid %d WAITS  (t+dt=%.8f > dependenciesMinTime=%.8f)\n",
               gid, endOfStepTime, getDependenciesMinTime());
#endif
        //mark this neuron as asleep waiting for a given min dependencies time
        dependenciesTimeNeuronWaitsFor = endOfStepTime;
        //release dependenciesLock and sleep until woken up by TimeDependencies::dependenciesWaitCondition
        libhpx_cond_wait(&this->dependenciesWaitCondition, &this->dependenciesLock);
    }
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
    printf("Neuron::waitForTimeDependencyNeurons: gid %d DONE   (t+dt=%.8f <= dependenciesMinTime=%.8f)\n",
           gid, endOfStepTime, getDependenciesMinTime());
#endif
    libhpx_mutex_unlock(&this->dependenciesLock);
}
