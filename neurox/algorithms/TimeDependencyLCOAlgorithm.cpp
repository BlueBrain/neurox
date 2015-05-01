#include "neurox/algorithms/TimeDependencyLCOAlgorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

floble_t DERIVED_CLASS_NAME::TimeDependencies::notificationIntervalRatio = 1;
double   DERIVED_CLASS_NAME::TimeDependencies::teps =1e-8;

DERIVED_CLASS_NAME::DERIVED_CLASS_NAME(){}

DERIVED_CLASS_NAME::~DERIVED_CLASS_NAME() {}

const AlgorithmType DERIVED_CLASS_NAME::getType()
{
    return AlgorithmType::BackwardEulerTimeDependencyLCO;
}

const char* DERIVED_CLASS_NAME::getTypeString()
{
    return "BackwardEulerTimeDependencyLCO";
}

void DERIVED_CLASS_NAME::Init()
{
    if (inputParams->allReduceAtLocality)
        throw std::runtime_error("Cant run BackwardEulerTimeDependencyLCO with allReduceAtLocality\n");

    const int allReducesCount = 0;
    hpx_bcast_rsync(AllReduceAlgorithm::AllReducesInfo::SetReductionsPerCommStep,
                    &allReducesCount, sizeof(int));
}

void DERIVED_CLASS_NAME::Clear() {}

double DERIVED_CLASS_NAME::Launch()
{
    int totalSteps = Algorithm::getTotalStepsCount();
    hpx_time_t now = hpx_time_now();
    neurox_hpx_call_neurons_lco(Branch::BackwardEuler, &totalSteps, sizeof(int));
    double elapsedTime = hpx_time_elapsed_ms(now)/1e3;
    input::Debugger::RunCoreneuronAndCompareAllBranches();
    return elapsedTime;
}

void DERIVED_CLASS_NAME::Run(Branch* b, const void* args)
{
    int steps = *(int*)args;

    if (b->soma)
    {
      TimeDependencies * timeDependencies = (TimeDependencies*) b->soma->algorithmMetaData;

      //fixes crash for Algorithm::All when TimeDependency algorithm starts at t=inputParams->tend*2
      //increase notification and dependencies time
      for (Neuron::Synapse *& s : b->soma->synapses)
          s->nextNotificationTime += b->nt->_t;
      timeDependencies->IncreseDependenciesTime(b->nt->_t);
    }

    for (int step=0; step<steps; step++)
        b->BackwardEulerStep();
    // Input::Coreneuron::Debugger::stepAfterStepBackwardEuler(local, &nrn_threads[this->nt->id], secondorder); //SMP ONLY

#ifndef NDEBUG
        printf("-- neuron %d finished\n", b->soma->gid);
#endif
}

void DERIVED_CLASS_NAME::StepBegin(Branch* b)
{
    if (b->soma)
    {
        TimeDependencies * timeDependencies = (TimeDependencies*) b->soma->algorithmMetaData;
        //inform time dependants that must be notified in this step
        timeDependencies->SendSteppingNotification(b->nt->_t, b->nt->_dt, b->soma->gid, b->soma->synapses);
        //wait until Im sure I can start and finalize this step at t+dt
        timeDependencies->WaitForTimeDependencyNeurons(b->nt->_t, b->nt->_dt, b->soma->gid);
    }
}

void DERIVED_CLASS_NAME::StepEnd(Branch* b, hpx_t)
{
    input::Debugger::SingleNeuronStepAndCompare(&nrn_threads[b->nt->id], b, inputParams->secondorder);
}

void DERIVED_CLASS_NAME::AfterReceiveSpikes(
        Branch *b, hpx_t target, neuron_id_t preNeuronId,
        spike_time_t spikeTime, spike_time_t maxTime)
{
    //inform soma of this neuron of new time dependency update
    hpx_t topBranchAddr = b->soma ? target : b->branchTree->topBranchAddr;
    if (b->soma)
    {
        TimeDependencies * timeDependencies = (TimeDependencies*) b->soma->algorithmMetaData;
        timeDependencies->UpdateTimeDependency(preNeuronId, maxTime);
    }
    else
        hpx_call(topBranchAddr, Branch::UpdateTimeDependency, HPX_NULL,
             &preNeuronId, sizeof(neuron_id_t), &maxTime, sizeof(spike_time_t));
}

hpx_t DERIVED_CLASS_NAME::SendSpikes(Neuron* neuron, double tt, double t)
{
    const floble_t notifRatio = TimeDependencyLCOAlgorithm::TimeDependencies::notificationIntervalRatio;
    const double teps = TimeDependencyLCOAlgorithm::TimeDependencies::teps;

    for (Neuron::Synapse *& s : neuron->synapses)
    {
        s->nextNotificationTime     = t+(s->minDelay+neuron->refractoryPeriod)*notifRatio;
        spike_time_t maxTimeAllowed = t+teps+s->minDelay+neuron->refractoryPeriod;

        hpx_lco_wait_reset(s->previousSpikeLco); //reset LCO to be used next
        //any spike or step notification happening after must wait for this spike delivery

        hpx_call(s->branchAddr, Branch::AddSpikeEvent, s->previousSpikeLco,
            &neuron->gid, sizeof(neuron_id_t), &tt, sizeof(spike_time_t),
            &maxTimeAllowed, sizeof(spike_time_t));

#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
        printf("Neuron::sendSpikes: gid %d at time %.3f, informs gid %d of next notif time =%.3f\n",
               this->gid, tt, s->destinationGid, t, s->nextNotificationTime);
#endif
    }
    return HPX_NULL;
}



DERIVED_CLASS_NAME::TimeDependencies::TimeDependencies()
{
    libhpx_cond_init(&this->dependenciesWaitCondition);
    libhpx_mutex_init(&this->dependenciesLock);
    this->dependenciesTimeNeuronWaitsFor = 0; //0 means not waiting
}

DERIVED_CLASS_NAME::TimeDependencies::~TimeDependencies()
{
    libhpx_cond_destroy(&this->dependenciesWaitCondition);
    libhpx_mutex_destroy(&this->dependenciesLock);
}

size_t DERIVED_CLASS_NAME::TimeDependencies::GetDependenciesCount()
{
    size_t size = -1;
    libhpx_mutex_lock(&this->dependenciesLock);
    size = dependenciesMap.size();
    libhpx_mutex_unlock(&this->dependenciesLock);
    return size;
}

void DERIVED_CLASS_NAME::TimeDependencies::IncreseDependenciesTime(floble_t t)
{
    libhpx_mutex_lock(&this->dependenciesLock);
    for (auto & dependency : dependenciesMap)
        dependency.second += t;
    libhpx_mutex_unlock(&this->dependenciesLock);
}

floble_t DERIVED_CLASS_NAME::TimeDependencies::GetDependenciesMinTime()
{
    assert(dependenciesMap.size()>0);
    return std::min_element(
        dependenciesMap.begin(), dependenciesMap.end(),
        [] (pair<neuron_id_t, floble_t> const& lhs, pair<neuron_id_t, floble_t> const& rhs)
            {return lhs.second < rhs.second;} )->second;
}

void DERIVED_CLASS_NAME::TimeDependencies::UpdateTimeDependency(
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
          if (GetDependenciesMinTime() >= dependenciesTimeNeuronWaitsFor) //and new min time allows neuron to proceed
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

void DERIVED_CLASS_NAME::TimeDependencies::WaitForTimeDependencyNeurons(floble_t t, floble_t dt, int gid)
{
    //if I have no dependencies... I'm free to go!
    if (dependenciesMap.size()==0) return;

#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
    printf("== %d enters TimeDependencies::waitForTimeDependencyNeurons\n", gid);
#endif
    libhpx_mutex_lock(&this->dependenciesLock);
    if (GetDependenciesMinTime() < t+dt) //if I cant proceed
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
    assert(GetDependenciesMinTime()>=t+dt);
    libhpx_mutex_unlock(&this->dependenciesLock);
#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
    printf("== %d leaves TimeDependencies::waitForTimeDependencyNeurons\n", gid);
#endif
}

void DERIVED_CLASS_NAME::TimeDependencies::SendSteppingNotification(floble_t t, floble_t dt, int gid, std::vector<Neuron::Synapse*> & synapses)
{
    for (Neuron::Synapse *& s : synapses)
        if (s->nextNotificationTime-teps <= t+dt) //if in this time step
            //(-teps to give or take few nanosecs for correction of floating point time roundings)
        {
            assert(s->nextNotificationTime >= t); //must have been covered by previous steps
            s->nextNotificationTime     = t+s->minDelay*TimeDependencies::notificationIntervalRatio;
            spike_time_t maxTimeAllowed = t+TimeDependencies::teps+s->minDelay;

            //wait for previous synapse to be delivered (if any) before telling post-syn neuron to proceed in time
            hpx_lco_wait(s->previousSpikeLco);
            hpx_call(s->topBranchAddr, Branch::UpdateTimeDependency, HPX_NULL,
                     &gid, sizeof(neuron_id_t), &maxTimeAllowed, sizeof(spike_time_t));

#if !defined(NDEBUG) && defined(PRINT_TIME_DEPENDENCY)
            printf("## %d notifies %d he can proceed up to %.6fms\n", gid, s->destinationGid, maxTimeAllowed);
#endif
        }
}
