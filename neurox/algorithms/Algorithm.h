#pragma once
#include "neurox.h"

namespace neurox
{

namespace algorithms
{

enum AlgorithmType {
    BackwardEulerDebugMode=0,
    BackwardEulerAllReduce=1,
    BackwardEulerSlidingTimeWindow=2,
    BackwardEulerTimeDependencyLCO=3,
    All=9 //All except debug (for benchmarking)
};

/**
 * @brief The Algorithm class
 * Virtual class describing all methods required to be implemented by all different algorithms
 */
class Algorithm
{
  public:
    Algorithm() {};
    virtual ~Algorithm() {};

    static Algorithm* New(AlgorithmType);
    virtual AlgorithmType getType() = 0;

    virtual void init(Branch*) {};
    virtual void clear(Branch*) {};
    virtual void stepBegin(Branch*) {};
    virtual void stepEnd(Branch*) {};
    virtual void commStepBegin(Branch*) {};
    virtual void commStepEnd(Branch*) {};
    virtual void afterSpike(Branch*) {};
};

}; //algorithms

}; //neurox

#include "neurox/algorithms/BackwardEulerAllReduceAlgorithm.h"
#include "neurox/algorithms/BackwardEulerDebugModeAlgorithm.h"
#include "neurox/algorithms/BackwardEulerSlidingTimeWindowAlgorithm.h"
#include "neurox/algorithms/BackwardEulerTimeDependencyLCOAlgorithm.h"
