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
