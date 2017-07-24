#pragma once
#include "neurox.h"

#define DERIVED_CLASS_NAME BackwardEulerSlidingTimeWindowAlgorithm

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

    const AlgorithmType getType() override;
    const char* getTypeString() override;

    void Init() override;
    void Clear() override;
    double Launch() override;

    void StepBegin(Branch*) override;
    void StepEnd(Branch*, hpx_t) override;
    void CommStepBegin(Branch*) override;
    void CommStepEnd(Branch*) override;
    hpx_t SendSpikes(Neuron*, double, double) override;

    const size_t allReducesCount = 2;
    static hpx_t * allReduces;

  private:
};

}; //algorithm

}; //neurox
