#include "neurox/algorithms/BackwardEulerSlidingTimeWindowAlgorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

hpx_t* DERIVED_CLASS_NAME::allReduces = nullptr;

DERIVED_CLASS_NAME::DERIVED_CLASS_NAME()
{
    Neuron::SlidingTimeWindow::reductionsPerCommStep = DERIVED_CLASS_NAME::allReducesCount;
}

DERIVED_CLASS_NAME::~DERIVED_CLASS_NAME() {}

const AlgorithmType DERIVED_CLASS_NAME::getType()
{
    return AlgorithmType::BackwardEulerSlidingTimeWindow;
}

const char* DERIVED_CLASS_NAME::getTypeString()
{
    return "BackwardEulerSlidingTimeWindow";
}

void DERIVED_CLASS_NAME::Init() {
    BackwardEulerAllReduceAlgorithm::SubscribeAllReduces(
                DERIVED_CLASS_NAME::allReduces,
                DERIVED_CLASS_NAME::allReducesCount);
}

void DERIVED_CLASS_NAME::Finalize() {
    BackwardEulerAllReduceAlgorithm::UnsubscribeAllReduces(
                DERIVED_CLASS_NAME::allReduces,
                DERIVED_CLASS_NAME::allReducesCount);
}

double DERIVED_CLASS_NAME::Run()
{
    int totalSteps = Algorithm::getTotalStepsCount();
    hpx_time_t now = hpx_time_now();
    if (inputParams->allReduceAtLocality)
        hpx_bcast_rsync(Branch::BackwardEulerOnLocality, &totalSteps, sizeof(int));
    else
        neurox_hpx_call_neurons_lco(Branch::BackwardEuler, &totalSteps, sizeof(int));
    double elapsedTime = hpx_time_elapsed_ms(now)/1e3;
    input::Debugger::RunCoreneuronAndCompareAllBranches();
    return elapsedTime;
}

void DERIVED_CLASS_NAME::StepBegin(Branch*) {}

void DERIVED_CLASS_NAME::StepEnd(Branch* b, hpx_t spikesLco)
{
    BackwardEulerAllReduceAlgorithm::WaitForSpikesDelivery(b, spikesLco);
    input::Debugger::SingleNeuronStepAndCompare(&nrn_threads[b->nt->id], b, inputParams->secondorder);
}

void DERIVED_CLASS_NAME::CommStepBegin(Branch*) {}

void DERIVED_CLASS_NAME::CommStepEnd(Branch*) {}

hpx_t DERIVED_CLASS_NAME::SendSpikes(Neuron* neuron, double tt, double)
{
    return BackwardEulerAllReduceAlgorithm::SendSpikes2(neuron,tt);
}
