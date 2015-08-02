#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::algorithms;

Algorithm* Algorithm::New(AlgorithmType type)
{
    switch (type)
    {
    case AlgorithmType::BackwardEulerCoreneuronDebug:
        return new CoreneuronDebugAlgorithm();
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
           getTypeString(), neurons->size(), inputParams->tstop/1000, inputParams->dt);
    fflush(stdout);
}

int Algorithm::getTotalStepsCount()
{
    return (inputParams->tstop - inputParams->tstart) / inputParams->dt;
}
