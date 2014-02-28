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

    const AlgorithmType getType() override;
    const char* getTypeString() override;

    void Init() override;
    void Finalize() override;
    double Run() override;

    void StepBegin(Branch*) override;
    void StepEnd(Branch*, hpx_t) override;
    void CommStepBegin(Branch*) override;
    void CommStepEnd(Branch*) override;
    hpx_t SendSpikes(Neuron* b, double tt, double t) override;
};

}; //algorithm

}; //neurox
