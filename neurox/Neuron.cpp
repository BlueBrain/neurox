#include "neurox/neurox.h"
#include <cstring>
#include <algorithm>
#include <utility>

using namespace neurox;
using namespace neurox::Solver;

Neuron::Neuron(neuron_id_t neuronId, floble_t APthreshold):
    gid(neuronId), threshold(APthreshold),
    timeDependencies(nullptr), commBarrier(nullptr), slidingTimeWindow(nullptr)
{
    this->synapsesTransmissionFlag = false;
    this->synapsesMutex = hpx_lco_sema_new(1);
    this->refractoryPeriod=0;
    this->commBarrier =
            inputParams->algorithm == Algorithm::BackwardEulerDebugWithCommBarrier
            ? new CommunicationBarrier() : NULL;
    this->timeDependencies =
            inputParams->algorithm == Algorithm::ALL ||
            inputParams->algorithm == Algorithm::BackwardEulerWithTimeDependencyLCO
            ? new TimeDependencies() : NULL;
    this->slidingTimeWindow =
            inputParams->algorithm == Algorithm::ALL ||
            inputParams->algorithm == Algorithm::BackwardEulerWithSlidingTimeWindow ||
            inputParams->algorithm == Algorithm::BackwardEulerWithAllReduceBarrier
            ? new SlidingTimeWindow() : NULL;
    assert(TimeDependencies::notificationIntervalRatio>0 && TimeDependencies::notificationIntervalRatio<=1);
    assert(Neuron::CommunicationBarrier::commStepSize % Neuron::SlidingTimeWindow::reductionsPerCommStep==0);
}

Neuron::~Neuron()
{
    if (synapsesMutex!=HPX_NULL)
        hpx_lco_delete_sync(synapsesMutex);
    for (Synapse *& s : synapses)
        delete s;
    delete timeDependencies;
    delete commBarrier;
    delete slidingTimeWindow;
}

Neuron::Synapse::Synapse(hpx_t branchAddr, floble_t minDelay, hpx_t topBranchAddr, int destinationGid)
    :branchAddr(branchAddr),minDelay(minDelay), topBranchAddr(topBranchAddr)
{
    this->nextNotificationTime=inputParams->tstart+TimeDependencies::teps+this->minDelay*TimeDependencies::notificationIntervalRatio;
    this->previousSpikeLco = hpx_lco_future_new(0);
    hpx_lco_set_rsync(this->previousSpikeLco, 0, NULL); //starts as set and will be reset when synapses happen
#ifndef NDEBUG
    this->destinationGid=destinationGid;
#endif
}

Neuron::Synapse::~Synapse()
{
    if (previousSpikeLco != HPX_NULL)
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

hpx_t Neuron::sendSpikes(floble_t t) //netcvode.cpp::PreSyn::send()
{
    if (synapses.size() == 0) return HPX_NULL; 

    spike_time_t tt = (spike_time_t) t+1e-10; //Coreneuron logic, do not change!
#if !defined(NDEBUG)
        printf("== Neuron gid %d spiked at %.3f ms\n", this->gid, tt);
#endif

    if (inputParams->algorithm==neurox::Algorithm::BackwardEulerDebugWithCommBarrier)
    {
        if (this->commBarrier->allSpikesLco == HPX_NULL) //first use
            this->commBarrier->allSpikesLco = hpx_lco_and_new(synapses.size());
        else
            hpx_lco_reset_sync(this->commBarrier->allSpikesLco); //reset to use after

        for (Synapse *& s : synapses)
            //deliveryTime (t+delay) is handled on post-syn side (diff value for every NetCon)
            hpx_call(s->branchAddr, Branch::addSpikeEvent, this->commBarrier->allSpikesLco,
                &this->gid, sizeof(neuron_id_t), &tt, sizeof(spike_time_t));
    }
    else if (inputParams->algorithm==neurox::Algorithm::BackwardEulerWithSlidingTimeWindow
          || inputParams->algorithm==neurox::Algorithm::BackwardEulerWithAllReduceBarrier)
    {
        hpx_t newSynapsesLco = hpx_lco_and_new(synapses.size());
        for (Synapse *& s : synapses)
            hpx_call(s->branchAddr, Branch::addSpikeEvent, newSynapsesLco,
                &this->gid, sizeof(neuron_id_t), &tt, sizeof(spike_time_t));
        return newSynapsesLco;
    }
    else if (inputParams->algorithm==neurox::Algorithm::BackwardEulerWithTimeDependencyLCO)
    {
        for (Synapse *& s : synapses)
        {
            s->nextNotificationTime     = t+(s->minDelay+refractoryPeriod)*TimeDependencies::notificationIntervalRatio;
            spike_time_t maxTimeAllowed = t+TimeDependencies::teps+s->minDelay+refractoryPeriod;

            hpx_lco_wait_reset(s->previousSpikeLco); //reset LCO to be used next
            //any spike or step notification happening after must wait for this spike delivery

            hpx_call(s->branchAddr, Branch::addSpikeEvent, s->previousSpikeLco,
                &this->gid, sizeof(neuron_id_t), &tt, sizeof(spike_time_t),
                &maxTimeAllowed, sizeof(spike_time_t));

#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
            printf("Neuron::sendSpikes: gid %d at time %.3f, informs gid %d of next notif time =%.3f\n",
                   this->gid, tt, s->destinationGid, t, s->nextNotificationTime);
#endif
        }
    }
    return HPX_NULL;
}

/////////////////// Neuron::CommunicationBarrier ///////////////////

Neuron::CommunicationBarrier::CommunicationBarrier()
{
    this->allSpikesLco = HPX_NULL;
    assert(Neuron::CommunicationBarrier::commStepSize <= 4);
}

Neuron::CommunicationBarrier::~CommunicationBarrier()
{
    if (allSpikesLco != HPX_NULL)
        hpx_lco_delete_sync(allSpikesLco);
}

//////////////////// Neuron::SlidingTimeWindow ///////////////////////

Neuron::SlidingTimeWindow::SlidingTimeWindow()
{
    for (int s=0; s< CommunicationBarrier::commStepSize-1; s++)
        this->spikesLcoQueue.push(HPX_NULL);
}

Neuron::SlidingTimeWindow::~SlidingTimeWindow()
{
    for (int i=0; i<spikesLcoQueue.size(); i++)
    {
        hpx_t queuedSpikesLco = spikesLcoQueue.front();
        if (queuedSpikesLco != HPX_NULL)
            hpx_lco_delete_sync(queuedSpikesLco);
        spikesLcoQueue.pop();
    }
}

hpx_action_t Neuron::SlidingTimeWindow::subscribeAllReduce = 0;
int Neuron::SlidingTimeWindow::subscribeAllReduce_handler(const hpx_t * allreduces, const size_t size)
{
    neurox_hpx_pin(Branch);
    SlidingTimeWindow * stw = local->soma->slidingTimeWindow;
    stw->allReduceFuture = new hpx_t[SlidingTimeWindow::reductionsPerCommStep];
    stw->allReduceLco = new hpx_t[SlidingTimeWindow::reductionsPerCommStep];
    stw->allReduceId = new int[SlidingTimeWindow::reductionsPerCommStep];
    for (int i=0; i<size/sizeof(hpx_t); i++)
    {
        stw->allReduceLco[i] = allreduces[i];
        stw->allReduceFuture[i] = hpx_lco_future_new(0); //no value to be reduced
        stw->allReduceId[i] = hpx_process_collective_allreduce_subscribe(
                allreduces[i], hpx_lco_set_action, stw->allReduceFuture[i]);
    }
    neurox_hpx_unpin;
}

hpx_action_t Neuron::SlidingTimeWindow::unsubscribeAllReduce = 0;
int Neuron::SlidingTimeWindow::unsubscribeAllReduce_handler(const hpx_t * allreduces, const size_t size)
{
    neurox_hpx_pin(Branch);
    SlidingTimeWindow * stw = local->soma->slidingTimeWindow;
    for (int i=0; i<size/sizeof(hpx_t); i++)
    {
      hpx_process_collective_allreduce_unsubscribe(allreduces[i], stw->allReduceId[i]);
      if (stw->allReduceFuture[i]!= HPX_NULL)
          hpx_lco_delete_sync(stw->allReduceFuture[i]);
    }
    delete [] stw->allReduceLco; stw->allReduceLco=nullptr;
    delete [] stw->allReduceFuture; stw->allReduceFuture=nullptr;
    delete [] stw->allReduceId; stw->allReduceId=nullptr;
    neurox_hpx_unpin;
}

int Neuron::SlidingTimeWindow::reductionsPerCommStep = -1;
std::vector<hpx_t>* Neuron::SlidingTimeWindow::AllReduceLocality::localityNeurons = nullptr;
hpx_t* Neuron::SlidingTimeWindow::AllReduceLocality::allReduceFuture = nullptr;
hpx_t* Neuron::SlidingTimeWindow::AllReduceLocality::allReduceLco = nullptr;
int* Neuron::SlidingTimeWindow::AllReduceLocality::allReduceId = nullptr;

hpx_action_t Neuron::SlidingTimeWindow::setReductionsPerCommStep = 0;
int Neuron::SlidingTimeWindow::setReductionsPerCommStep_handler(const int* val, const size_t)
{
    neurox_hpx_pin(uint64_t);
    reductionsPerCommStep = *val;
    neurox_hpx_unpin;
}

hpx_action_t Neuron::SlidingTimeWindow::AllReduceLocality::subscribeAllReduce = 0;
int Neuron::SlidingTimeWindow::AllReduceLocality::subscribeAllReduce_handler(const hpx_t * allreduces, const size_t size)
{
    neurox_hpx_pin(uint64_t);
    assert(inputParams->allReduceAtLocality);
    AllReduceLocality::allReduceLco = new hpx_t[SlidingTimeWindow::reductionsPerCommStep];
    AllReduceLocality::allReduceFuture = new hpx_t[SlidingTimeWindow::reductionsPerCommStep];
    AllReduceLocality::allReduceId = new int[SlidingTimeWindow::reductionsPerCommStep];
    for (int i=0; i<size/sizeof(hpx_t); i++)
    {
        allReduceLco[i] = allreduces[i];
        allReduceFuture[i] = hpx_lco_future_new(0); //no value to be reduced
        allReduceId[i] = hpx_process_collective_allreduce_subscribe(
                allreduces[i], hpx_lco_set_action, allReduceFuture[i]);
    }
    neurox_hpx_unpin;
}

hpx_action_t Neuron::SlidingTimeWindow::AllReduceLocality::unsubscribeAllReduce = 0;
int Neuron::SlidingTimeWindow::AllReduceLocality::unsubscribeAllReduce_handler(const hpx_t * allreduces, const size_t size)
{
    neurox_hpx_pin(uint64_t);
    assert(inputParams->allReduceAtLocality);
    for (int i=0; i<size/sizeof(hpx_t); i++)
    {
      hpx_process_collective_allreduce_unsubscribe(allreduces[i], allReduceId[i]);
      hpx_lco_delete_sync(allReduceFuture[i]);
    }
    delete [] allReduceLco; allReduceLco=nullptr;
    delete [] allReduceFuture; allReduceFuture=nullptr;
    delete [] allReduceId; allReduceId=nullptr;
    neurox_hpx_unpin;
}

hpx_action_t Neuron::SlidingTimeWindow::init = 0;
void Neuron::SlidingTimeWindow::init_handler
    (const void*, const size_t) {}

hpx_action_t Neuron::SlidingTimeWindow::reduce = 0;
void Neuron::SlidingTimeWindow::reduce_handler
    (void* , const void* , const size_t) {}

////////////////////// Neuron::TimeDependencies ///////////////////

Neuron::TimeDependencies::TimeDependencies()
{
    libhpx_cond_init(&this->dependenciesWaitCondition);
    libhpx_mutex_init(&this->dependenciesLock);
    this->dependenciesTimeNeuronWaitsFor = 0; //0 means not waiting
}

Neuron::TimeDependencies::~TimeDependencies()
{
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

void Neuron::TimeDependencies::increseDependenciesTime(floble_t t)
{
    libhpx_mutex_lock(&this->dependenciesLock);
    for (auto & dependency : dependenciesMap)
        dependency.second += t;
    libhpx_mutex_unlock(&this->dependenciesLock);
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
        neuron_id_t srcGid, floble_t dependencyNotificationTime, neuron_id_t myGid,  bool initializationPhase)
{
    libhpx_mutex_lock(&this->dependenciesLock);
    if (initializationPhase)
    {
        if (dependenciesMap.find(srcGid) == dependenciesMap.end()) //in execution without branching this is always true
          dependenciesMap[srcGid] = dependencyNotificationTime;
        else
          dependenciesMap.at(srcGid) = std::max(dependenciesMap.at(srcGid), dependencyNotificationTime);
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
            spike_time_t maxTimeAllowed = t+TimeDependencies::teps+s->minDelay;

            //wait for previous synapse to be delivered (if any) before telling post-syn neuron to proceed in time
            hpx_lco_wait(s->previousSpikeLco);
            hpx_call(s->topBranchAddr, Branch::updateTimeDependency, HPX_NULL,
                     &gid, sizeof(neuron_id_t), &maxTimeAllowed, sizeof(spike_time_t));

#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
            printf("## %d notifies %d he can proceed up to %.6fms\n", gid, s->destinationGid, maxTimeAllowed);
#endif
        }
}

void Neuron::registerHpxActions()
{
    neurox_hpx_register_action(neurox_single_var_action, SlidingTimeWindow::subscribeAllReduce);
    neurox_hpx_register_action(neurox_single_var_action, SlidingTimeWindow::unsubscribeAllReduce);
    neurox_hpx_register_action(neurox_single_var_action, SlidingTimeWindow::AllReduceLocality::subscribeAllReduce);
    neurox_hpx_register_action(neurox_single_var_action, SlidingTimeWindow::AllReduceLocality::unsubscribeAllReduce);
    neurox_hpx_register_action(neurox_single_var_action, SlidingTimeWindow::setReductionsPerCommStep);
    neurox_hpx_register_action(neurox_reduce_op_action,  SlidingTimeWindow::init);
    neurox_hpx_register_action(neurox_reduce_op_action,  SlidingTimeWindow::reduce);
}
