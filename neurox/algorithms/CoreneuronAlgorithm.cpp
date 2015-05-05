#include "neurox/algorithms/CoreneuronAlgorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

DERIVED_CLASS_NAME::DERIVED_CLASS_NAME() {}

DERIVED_CLASS_NAME::~DERIVED_CLASS_NAME() {}

DERIVED_CLASS_NAME::CommunicationBarrier::CommunicationBarrier()
{
    assert(0);
}

DERIVED_CLASS_NAME::CommunicationBarrier::~CommunicationBarrier()
{
    assert(0);
}

const AlgorithmType DERIVED_CLASS_NAME::getType()
{
    return AlgorithmType::BackwardEulerCoreneuron;
}

const char* DERIVED_CLASS_NAME::getTypeString()
{
    return "BackwardEulerCoreneuron";
}

void DERIVED_CLASS_NAME::Init()
{
    assert(0);
}

void DERIVED_CLASS_NAME::Clear() {}

double DERIVED_CLASS_NAME::Launch()
{
    int commStepSize = CoreneuronAlgorithm::CommunicationBarrier::commStepSize;
    int totalSteps = Algorithm::getTotalStepsCount();
    hpx_time_t now = hpx_time_now();
    assert(0);
    double elapsedTime = hpx_time_elapsed_ms(now)/1e3;
    input::Debugger::CompareAllBranches();
    return elapsedTime;
}

void DERIVED_CLASS_NAME::StepBegin(Branch*)
{
    assert(0);
}

void DERIVED_CLASS_NAME::StepEnd(Branch* b, hpx_t)
{
    assert(0);
}

void DERIVED_CLASS_NAME::Run(Branch* b, const void* args)
{
    assert(0);
}

hpx_t DERIVED_CLASS_NAME::SendSpikes(Neuron* neuron, double tt, double)
{
    assert(0);
    return HPX_NULL;
}
