#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::algorithms;

Algorithm* Algorithm::New(AlgorithmType type)
{
    switch (type)
    {
    case AlgorithmType::BackwardEulerDebugMode:
        return new BackwardEulerDebugModeAlgorithm();
    case AlgorithmType::BackwardEulerAllReduce:
        return new BackwardEulerAllReduceAlgorithm();
    case AlgorithmType::BackwardEulerSlidingTimeWindow:
        return new BackwardEulerSlidingTimeWindowAlgorithm();
    case AlgorithmType::BackwardEulerTimeDependencyLCO:
        return new BackwardEulerTimeDependencyLCOAlgorithm();
    default:
        return nullptr;
    }
    return nullptr;
}

void Algorithm::PrintStartInfo()
{
    printf("neurox::Algorithm::%s (%d neurons, t=%.03f secs, dt=%.03f milisecs\n",
           getTypeString(), neurons->size(), inputParams->tstop/1000, inputParams->dt);
    fflush(stdout);
}

int Algorithm::getTotalStepsCount()
{
    return (inputParams->tstop - inputParams->tstart) / inputParams->dt;
}
