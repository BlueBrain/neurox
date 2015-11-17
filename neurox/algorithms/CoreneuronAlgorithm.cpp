#include "neurox/algorithms/CoreneuronAlgorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

CoreneuronAlgorithm::CoreneuronAlgorithm() {}

CoreneuronAlgorithm::~CoreneuronAlgorithm() {}

CoreneuronAlgorithm::CommunicationBarrier::CommunicationBarrier()
{
    assert(0);
}

CoreneuronAlgorithm::CommunicationBarrier::~CommunicationBarrier()
{
    assert(0);
}

const AlgorithmType CoreneuronAlgorithm::getType()
{
    return AlgorithmType::BackwardEulerCoreneuron;
}

const char* CoreneuronAlgorithm::getTypeString()
{
    return "BackwardEulerCoreneuron";
}

void CoreneuronAlgorithm::Init()
{
    assert(0);
}

void CoreneuronAlgorithm::Clear() {}

double CoreneuronAlgorithm::Launch()
{
    int commStepSize = CoreneuronAlgorithm::CommunicationBarrier::commStepSize;
    int totalSteps = Algorithm::getTotalStepsCount();
    hpx_time_t now = hpx_time_now();
    assert(0);
    double elapsedTime = hpx_time_elapsed_ms(now)/1e3;
    input::Debugger::CompareAllBranches();
    return elapsedTime;
}

void CoreneuronAlgorithm::StepBegin(Branch*)
{
    assert(0);
}

void CoreneuronAlgorithm::StepEnd(Branch* b, hpx_t)
{
    assert(0);
}

void CoreneuronAlgorithm::Run(Branch* b, const void* args)
{
    assert(0);
}

hpx_t CoreneuronAlgorithm::SendSpikes(Neuron* neuron, double tt, double)
{
    assert(0);
    return HPX_NULL;
}
