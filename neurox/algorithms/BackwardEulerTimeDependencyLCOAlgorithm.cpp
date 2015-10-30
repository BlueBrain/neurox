#include "neurox/algorithms/BackwardEulerTimeDependencyLCOAlgorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

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
    hpx_bcast_rsync(neurox::Neuron::SlidingTimeWindow::SetReductionsPerCommStep,
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

void DERIVED_CLASS_NAME::Run(Branch*) {}

void DERIVED_CLASS_NAME::StepBegin(Branch* b)
{
    if (b->soma)
    {
        //inform time dependants that must be notified in this step
        b->soma->timeDependencies->SendSteppingNotification(b->nt->_t, b->nt->_dt, b->soma->gid, b->soma->synapses);
        //wait until Im sure I can start and finalize this step at t+dt
        b->soma->timeDependencies->WaitForTimeDependencyNeurons(b->nt->_t, b->nt->_dt, b->soma->gid);
    }
}

void DERIVED_CLASS_NAME::StepEnd(Branch* b, hpx_t)
{
    input::Debugger::SingleNeuronStepAndCompare(&nrn_threads[b->nt->id], b, inputParams->secondorder);
}

void DERIVED_CLASS_NAME::AfterReceiveSpikes(
        Branch *local, hpx_t target, neuron_id_t preNeuronId,
        spike_time_t spikeTime, spike_time_t maxTime)
{
    //inform soma of this neuron of new time dependency update
    hpx_t topBranchAddr = local->soma ? target : local->branchTree->topBranchAddr;
    if (local->soma)
        local->soma->timeDependencies->UpdateTimeDependency(preNeuronId, maxTime);
    else
        hpx_call(topBranchAddr, Branch::UpdateTimeDependency, HPX_NULL,
             &preNeuronId, sizeof(neuron_id_t), &maxTime, sizeof(spike_time_t));
}

hpx_t DERIVED_CLASS_NAME::SendSpikes(Neuron* neuron, double tt, double t)
{
    for (Neuron::Synapse *& s : neuron->synapses)
    {
        s->nextNotificationTime     = t+(s->minDelay+neuron->refractoryPeriod)*Neuron::TimeDependencies::notificationIntervalRatio;
        spike_time_t maxTimeAllowed = t+Neuron::TimeDependencies::teps+s->minDelay+neuron->refractoryPeriod;

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
