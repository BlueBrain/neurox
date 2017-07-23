#pragma once
#include "neurox.h"

using namespace neurox;

namespace neurox
{

namespace algorithms
{

enum AlgorithmType {
    BackwardEulerDebugMode=0,
    BackwardEulerAllReduce=1,
    BackwardEulerSlidingTimeWindow=2,
    BackwardEulerTimeDependencyLCO=3,
    ALL=9 //All, except debug
};

class Algorithm
{
    constexpr virtual AlgorithmType getType() = 0;

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
