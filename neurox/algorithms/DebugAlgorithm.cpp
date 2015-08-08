#include "neurox/algorithms/DebugAlgorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

DERIVED_CLASS_NAME::DERIVED_CLASS_NAME() {}

DERIVED_CLASS_NAME::~DERIVED_CLASS_NAME() {}


DERIVED_CLASS_NAME::CommunicationBarrier::CommunicationBarrier()
{
    this->allSpikesLco = HPX_NULL;
    assert(CoreneuronAlgorithm::CommunicationBarrier::commStepSize <= 4);
}

DERIVED_CLASS_NAME::CommunicationBarrier::~CommunicationBarrier()
{
    if (allSpikesLco != HPX_NULL)
        hpx_lco_delete_sync(allSpikesLco);
}

const AlgorithmType DERIVED_CLASS_NAME::getType()
{
    return AlgorithmType::BackwardEulerDebug;
}

const char* DERIVED_CLASS_NAME::getTypeString()
{
    return "BackwardEulerCoreneuronDebug";
}

void DERIVED_CLASS_NAME::Init()
{
    const int allReducesCount = 0;
    hpx_bcast_rsync(AllReduceAlgorithm::AllReducesInfo::SetReductionsPerCommStep,
                    &allReducesCount, sizeof(int));
}

void DERIVED_CLASS_NAME::Clear() {}

double DERIVED_CLASS_NAME::Launch()
{
    int commStepSize = CoreneuronAlgorithm::CommunicationBarrier::commStepSize;
    int totalSteps = Algorithm::getTotalStepsCount();

    hpx_time_t now = hpx_time_now();
    for (int s=0; s<totalSteps; s+=commStepSize)
    {
        #ifdef NEUROX_TIME_STEPPING_VERBOSE
          if (hpx_get_my_rank()==0)
            DebugMessage(std::string("-- t="+std::to_string(inputParams->dt*s)+" ms\n").c_str());
        #endif

        //Reduction at locality is not implemented (this mode is for debugging only)
        neurox_hpx_call_neurons_lco(Branch::BackwardEuler, &commStepSize, sizeof(int));

        #ifndef NDEBUG
          if (neurox::ParallelExecution()) //if parallel execution... spike exchange
            hpx_bcast_rsync(neurox::input::Debugger::NrnSpikeExchange);
        #endif
    }
    double elapsedTime = hpx_time_elapsed_ms(now)/1e3;
    input::Debugger::CompareAllBranches();
    return elapsedTime;
}

void DERIVED_CLASS_NAME::StepBegin(Branch*) {}

void DERIVED_CLASS_NAME::StepEnd(Branch* b, hpx_t)
{
    input::Debugger::SingleNeuronStepAndCompare(&nrn_threads[b->nt->id], b, inputParams->secondorder);
}

void DERIVED_CLASS_NAME::Run(Branch* b, const void* args)
{
    int steps = *(int*)args;
    for (int step=0; step<steps; step++)
        b->BackwardEulerStep();
    // Input::Coreneuron::Debugger::stepAfterStepBackwardEuler(local, &nrn_threads[this->nt->id], secondorder); //SMP ONLY

    if (b->soma) //end of comm-step (steps is the number of steps per commSize)
    {
        CommunicationBarrier * commBarrier = (CommunicationBarrier *) b->soma->algorithmMetaData;
        if (commBarrier->allSpikesLco != HPX_NULL) //was set/used once
            hpx_lco_wait(commBarrier->allSpikesLco); //wait if needed
    }
}

hpx_t DERIVED_CLASS_NAME::SendSpikes(Neuron* neuron, double tt, double)
{
    CommunicationBarrier * commBarrier = (CommunicationBarrier *) neuron->algorithmMetaData;
    if (commBarrier->allSpikesLco == HPX_NULL) //first use
        commBarrier->allSpikesLco = hpx_lco_and_new(neuron->synapses.size());
    else
        hpx_lco_reset_sync(commBarrier->allSpikesLco); //reset to use after

    for (Neuron::Synapse *& s : neuron->synapses)
        //deliveryTime (t+delay) is handled on post-syn side (diff value for every NetCon)
        hpx_call(s->branchAddr, Branch::AddSpikeEvent, commBarrier->allSpikesLco,
            &neuron->gid, sizeof(neuron_id_t), &tt, sizeof(spike_time_t));

    return HPX_NULL;
}
