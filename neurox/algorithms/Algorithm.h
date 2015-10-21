#pragma once
#include "neurox.h"

using namespace  neurox;

namespace neurox
{

namespace algorithms
{

enum AlgorithmType {
    BackwardEulerDebugWithCommBarrier=0,
    BackwardEulerWithAllReduceBarrier=1,
    BackwardEulerWithSlidingTimeWindow=2,
    BackwardEulerWithTimeDependencyLCO=3,
    ALL=9 //All, except debug
};

class Algorithm
{
    virtual void init(Branch*) = 0;
    virtual void clear(Branch*) = 0;
    virtual void stepBegin(Branch*) = 0;
    virtual void stepEnd(Branch*) = 0;
    virtual void commStepBegin(Branch*) = 0;
    virtual void commStepEnd(Branch*) = 0;
    virtual void afterSpike(Branch*) = 0;
};

}; //algorithm

}; //neurox
