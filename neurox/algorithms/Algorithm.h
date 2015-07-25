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

    virtual void Run(AlgorithmType);

    /// Returns an instantiated class of the given type
    static Algorithm* New(AlgorithmType);

    const virtual AlgorithmType getType() = 0; ///> Returns class type
    const virtual char* getTypeString() = 0; ///> Returns class type string

    virtual void Init() {};
    virtual void Clear() {};
    virtual void StepBegin(Branch*) {};
    virtual void StepEnd(Branch*) {};
    virtual void CommStepBegin(Branch*) {};
    virtual void CommStepEnd(Branch*) {};
    virtual void AfterSpike(Branch*) {};
};

}; //algorithms

}; //neurox

#include "neurox/algorithms/BackwardEulerAllReduceAlgorithm.h"
#include "neurox/algorithms/BackwardEulerDebugModeAlgorithm.h"
#include "neurox/algorithms/BackwardEulerSlidingTimeWindowAlgorithm.h"
#include "neurox/algorithms/BackwardEulerTimeDependencyLCOAlgorithm.h"
