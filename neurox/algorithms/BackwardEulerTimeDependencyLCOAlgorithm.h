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
    void Init() override;
    void Clear() override;
    void StepBegin(Branch*) override;
    void StepEnd(Branch*) override;
    void CommStepBegin(Branch*) override;
    void CommStepEnd(Branch*) override;
    void AfterSpike(Branch*) override;
};

}; //algorithm

}; //neurox
