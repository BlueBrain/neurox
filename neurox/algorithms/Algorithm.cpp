#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::algorithms;

Algorithm* Algorithm::New(AlgorithmType type)
{
    switch (type)
    {
    case AlgorithmType::BackwardEulerDebugMode:
        return new BackwardEulerAllReduce();
        break;
    case AlgorithmType::BackwardEulerAllReduce:
        return new BackwardEulerAllReduce();
        break;
    case AlgorithmType::BackwardEulerSlidingTimeWindow:
        return new BackwardEulerSlidingTimeWindow();
        break;
    case AlgorithmType::BackwardEulerTimeDependencyLCO:
        return new BackwardEulerTimeDependencyLCO();
        break;
    default:
        return nullptr;
    }
}
