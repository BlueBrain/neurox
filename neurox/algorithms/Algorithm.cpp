#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::algorithms;

Algorithm* Algorithm::New(AlgorithmType type)
{
    switch (type)
    {
    case AlgorithmType::BackwardEulerDebug:
        return new DebugAlgorithm();
    case AlgorithmType::BackwardEulerCoreneuron:
        return new CoreneuronAlgorithm();
    case AlgorithmType::BackwardEulerAllReduce:
        return new AllReduceAlgorithm();
    case AlgorithmType::BackwardEulerSlidingTimeWindow:
        return new SlidingTimeWindowAlgorithm();
    case AlgorithmType::BackwardEulerTimeDependencyLCO:
        return new TimeDependencyLCOAlgorithm();
    default:
        return nullptr;
    }
    return nullptr;
};

void Algorithm::PrintStartInfo()
{
    printf("neurox::Algorithm::%s (%d neurons, t=%.03f secs, dt=%.03f milisecs\n",
           getTypeString(), neurox::neurons_count, input_params->tstop/1000, input_params->dt);
    fflush(stdout);
}

int Algorithm::getTotalStepsCount()
{
    return (input_params->tstop - input_params->tstart) / input_params->dt;
}
