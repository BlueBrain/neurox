#pragma once
#include "neurox.h"

#define DERIVED_CLASS_NAME BackwardEulerAllReduceAlgorithm

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

    static void SubscribeAllReduces  (hpx_t *& allReduces, size_t allReducesCount);
    static void UnsubscribeAllReduces(hpx_t *& allReduces, size_t allReducesCount);

    const size_t allReducesCount = 1;
    static hpx_t* allReduces;

  private:
};

}; //algorithm

}; //neurox
