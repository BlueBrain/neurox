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
    All=9, //All except debug (for benchmarking)
    None=-1 //If not instantiated
};

/**
 * @brief The Algorithm class
 * Virtual class describing all methods required to be implemented by all different algorithms
 */
class Algorithm
{
  public:
    Algorithm() = delete;
    static Algorithm * New(AlgorithmType);

    virtual AlgorithmType getType()
    {
        return AlgorithmType::None;
    };

    virtual void init(Branch*) = 0;
    virtual void clear(Branch*) = 0;
    virtual void stepBegin(Branch*) = 0;
    virtual void stepEnd(Branch*) = 0;
    virtual void commStepBegin(Branch*) = 0;
    virtual void commStepEnd(Branch*) = 0;
    virtual void afterSpike(Branch*) = 0;
};

}; //algorithms

}; //neurox

#include "neurox/algorithms/BackwardEulerAllReduceAlgorithm.h"
#include "neurox/algorithms/BackwardEulerDebugModeAlgorithm.h"
#include "neurox/algorithms/BackwardEulerSlidingTimeWindowAlgorithm.h"
#include "neurox/algorithms/BackwardEulerTimeDependencyLCOAlgorithm.h"
