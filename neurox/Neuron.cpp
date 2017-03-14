#include "neurox/neurox.h"
#include <cstring>
#include <algorithm>
#include <utility>

using namespace neurox;
using namespace neurox::Solver;

const double Neuron::teps = 1e-8;
const floble_t Neuron::TimeDependencies::notificationIntervalRatio = 1; //ratio (0,1]
const int Neuron::commStepSize = 4;

Neuron::Neuron(neuron_id_t neuronId, floble_t APthreshold, int thvar_index):
    gid(neuronId), threshold(APthreshold), thvar_index(thvar_index)
{
    this->synapsesTransmissionFlag = false;
    this->synapsesMutex = hpx_lco_sema_new(1);
    this->commStepAllSpikesLco = HPX_NULL;
    this->refractoryPeriod=0;
    this->timeDependencies = inputParams->algorithm == Algorithm::BackwardEulerWithPairwiseSteping ? new TimeDependencies() : NULL;
    //to avoid rounding errors, we send the notifications slighly before (so this should never be 1)
    assert(TimeDependencies::notificationIntervalRatio>0 && TimeDependencies::notificationIntervalRatio<=1);
    assert(Neuron::commStepSize <= 4);
}

Neuron::~Neuron()
{
    hpx_lco_delete_sync(synapsesMutex);
    for (Synapse *& s : synapses)
        delete s;
    delete timeDependencies;
    hpx_lco_delete_sync(commStepAllSpikesLco);
}

Neuron::Synapse::Synapse(hpx_t addr, floble_t minDelay, int destinationGid)
    :addr(addr),minDelay(minDelay)
{
    this->nextNotificationTime=inputParams->tstart+Neuron::teps+this->minDelay*TimeDependencies::notificationIntervalRatio;
    this->previousSpikeLco = hpx_lco_future_new(0);
    hpx_lco_set_rsync(this->previousSpikeLco, 0, NULL); //starts as set and will be reset when synapses happen
#ifndef NDEBUG
    this->destinationGid=destinationGid;
#endif
}

Neuron::Synapse::~Synapse()
{
    hpx_lco_delete_sync(previousSpikeLco);
}

size_t Neuron::getSynapseCount()
{
    size_t size = -1;
    hpx_lco_sema_p(synapsesMutex);
    size = synapses.size();
    hpx_lco_sema_v_sync(synapsesMutex);
    return size;
}

void Neuron::addSynapse(Synapse * syn)
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

void Neuron::sendSpikes(floble_t t, floble_t dt)
{
    if (synapses.size() == 0) return;
    spike_time_t tt = (spike_time_t) t+1e-10;

    //netcvode.cpp::PreSyn::send()
    if (inputParams->algorithm==neurox::Algorithm::BackwardEulerSyncFixedCommStepDebug
     || inputParams->algorithm==neurox::Algorithm::BackwardEulerAsyncFixedCommStep)
    {
        if (commStepAllSpikesLco == HPX_NULL) //first use
            commStepAllSpikesLco = hpx_lco_and_new(synapses.size());
        else
            hpx_lco_reset_sync(commStepAllSpikesLco); //reset to use after

        for (Synapse *& s : synapses)
            //deliveryTime (t+delay) is handled on post-syn side (diff value for every NetCon)
            hpx_call(s->addr, Branch::addSpikeEvent, commStepAllSpikesLco,
                &this->gid, sizeof(neuron_id_t), &tt, sizeof(spike_time_t));

        if (inputParams->algorithm==neurox::Algorithm::BackwardEulerSyncFixedCommStepDebug)
            hpx_lco_wait_reset(commStepAllSpikesLco); //Debug mode waits now, otherwise at the end of comm step
    }
    else if (inputParams->algorithm==neurox::Algorithm::BackwardEulerWithPairwiseSteping)
    {
        for (Synapse *& s : synapses)
        {
            s->nextNotificationTime     = t+(s->minDelay+refractoryPeriod)*TimeDependencies::notificationIntervalRatio;
            spike_time_t maxTimeAllowed = t+Neuron::teps+s->minDelay+refractoryPeriod;

            hpx_lco_wait_reset(s->previousSpikeLco); //reset LCO to be used next
            //any spike or step notification happening after must wait for this spike delivery

            hpx_call(s->addr, Branch::addSpikeEvent, s->previousSpikeLco,
                &this->gid, sizeof(neuron_id_t), &tt, sizeof(spike_time_t),
                &maxTimeAllowed, sizeof(spike_time_t));

#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
            printf("Neuron::sendSpikes: gid %d at time %.3f, informs gid %d of next notif time =%.3f\n",
                   this->gid, tt, s->destinationGid, t, s->nextNotificationTime);
#endif
        }
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
        neuron_id_t srcGid, neuron_id_t myGid, floble_t dependencyNotificationTime, bool initializationPhase)
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
        {
            dependenciesMap.at(srcGid) = dependencyNotificationTime;
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
            printf("-- %d (msg from %d) updates dependenciesMap(%d)=%.11f, notif time=%.11f, getDependenciesMinTime()=%.11f\n",
                   myGid, srcGid, srcGid, dependenciesMap.at(srcGid), dependencyNotificationTime, getDependenciesMinTime());
#endif
        }else
        {
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
            printf("-- %d (msg from %d) DOES NOT UPDATE dependenciesMap(%d)=%.11f, notif time=%.11f, getDependenciesMinTime()=%.11f\n",
                   myGid, srcGid, srcGid, dependenciesMap.at(srcGid), dependencyNotificationTime, getDependenciesMinTime());
#endif
        }

        if (dependenciesTimeNeuronWaitsFor > 0) //if neuron is waiting for a dependencies time update
          if (getDependenciesMinTime() >= dependenciesTimeNeuronWaitsFor) //and new min time allows neuron to proceed
          {
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
            printf("-- %d (msg from %d) wakes up producer, getDependenciesMinTime()=%.11f >= t+dt=%.11f\n",
                myGid, srcGid, getDependenciesMinTime(), dependenciesTimeNeuronWaitsFor);
#endif
            dependenciesTimeNeuronWaitsFor = 0; //mark neuron as not asleep anymore
            libhpx_cond_broadcast(&this->dependenciesWaitCondition); //wake up neuron
          }
    }
    libhpx_mutex_unlock(&this->dependenciesLock);
}

void Neuron::TimeDependencies::waitForTimeDependencyNeurons(floble_t t, floble_t dt, int gid)
{
    //if I have no dependencies... I'm free to go!
    if (dependenciesMap.size()==0) return;

#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
    printf("== %d enters TimeDependencies::waitForTimeDependencyNeurons\n", gid);
#endif
    libhpx_mutex_lock(&this->dependenciesLock);
    if (getDependenciesMinTime() < t+dt) //if I cant proceed
    {
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
        printf("== %d cant proceed and sleeps: getDependenciesMinTime()=%.11f < t+dt=%.11f\n", gid, getDependenciesMinTime(), t+dt);
#endif
        //mark this neuron as asleep waiting for a given min dependencies time
        dependenciesTimeNeuronWaitsFor = t+dt;
        //release dependenciesLock and sleep until woken up by TimeDependencies::dependenciesWaitCondition
        libhpx_cond_wait(&this->dependenciesWaitCondition, &this->dependenciesLock);
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
        printf("== %d wakes up: getDependenciesMinTime()=%.11f\n", gid, getDependenciesMinTime());
#endif
    }
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
    else
        printf("== %d proceeds: getDependenciesMinTime()=%.11f >= t+dt=%.11f\n", gid, getDependenciesMinTime(), t+dt);
#endif
    assert(getDependenciesMinTime()>=t+dt);
    libhpx_mutex_unlock(&this->dependenciesLock);
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
    printf("== %d leaves TimeDependencies::waitForTimeDependencyNeurons\n", gid);
#endif
}

void Neuron::TimeDependencies::sendSteppingNotification(floble_t t, floble_t dt, int gid, std::vector<Synapse*> & synapses)
{
    for (Synapse *& s : synapses)
        if (s->nextNotificationTime-teps <= t+dt) //if in this time step
            //(-teps to give or take few nanosecs for correction of floating point time roundings)
        {
            assert(s->nextNotificationTime >= t); //must have been covered by previous steps
            s->nextNotificationTime     = t+s->minDelay*TimeDependencies::notificationIntervalRatio;
            spike_time_t maxTimeAllowed = t+Neuron::teps+s->minDelay;

            //wait for previous synapse to be delivered (if any) before telling post-syn to proceed in time
            hpx_lco_wait(s->previousSpikeLco);
            hpx_call(s->addr, Branch::updateTimeDependencyValue, HPX_NULL,
                     &gid, sizeof(neuron_id_t), &maxTimeAllowed, sizeof(spike_time_t));

#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
            printf("## %d notifies %d he can proceed up to %.6fms\n", gid, s->destinationGid, maxTimeAllowed);
#endif
        }
}
