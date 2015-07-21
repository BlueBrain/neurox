#pragma once
#include "neurox.h"

using namespace neurox;

namespace neurox
{

namespace algorithms
{

class BackwardEulerTimeDependencyLCO : Algorithm
{
    constexpr virtual AlgorithmType getType() override
    {
        return AlgorithmType::BackwardEulerTimeDependencyLCO;
    };

    virtual void init(Branch*) override;
    virtual void clear(Branch*) override;
    virtual void stepBegin(Branch*) override;
    virtual void stepEnd(Branch*) override;
    virtual void commStepBegin(Branch*) override;
    virtual void commStepEnd(Branch*) override;
    virtual void afterSpike(Branch*) override;
};

}; //algorithm

}; //neurox
