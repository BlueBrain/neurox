#pragma once
#include "neurox.h"

#define DERIVED_CLASS_NAME BackwardEulerTimeDependencyLCOAlgorithm

using namespace neurox;

namespace neurox
{

namespace algorithms
{

class DERIVED_CLASS_NAME : public Algorithm
{
  public:
    DERIVED_CLASS_NAME();
    ~DERIVED_CLASS_NAME();

    AlgorithmType getType() override;
    void init(Branch*) override;
    void clear(Branch*) override;
    void stepBegin(Branch*) override;
    void stepEnd(Branch*) override;
    void commStepBegin(Branch*) override;
    void commStepEnd(Branch*) override;
    void afterSpike(Branch*) override;
};

}; //algorithm

}; //neurox
