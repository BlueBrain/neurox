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

void DERIVED_CLASS_NAME::Finalize() {}

double DERIVED_CLASS_NAME::Run()
{
    int totalSteps = Algorithm::getTotalStepsCount();
    hpx_time_t now = hpx_time_now();
    neurox_hpx_call_neurons_lco(Branch::BackwardEuler, &totalSteps, sizeof(int));
    double elapsedTime = hpx_time_elapsed_ms(now)/1e3;
    input::Debugger::RunCoreneuronAndCompareAllBranches();
    return elapsedTime;
}

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

void DERIVED_CLASS_NAME::StepEnd(Branch*, hpx_t) {}

void DERIVED_CLASS_NAME::CommStepBegin(Branch*) {}

void DERIVED_CLASS_NAME::CommStepEnd(Branch*) {}

void DERIVED_CLASS_NAME::AfterSpike(Branch*) {}
