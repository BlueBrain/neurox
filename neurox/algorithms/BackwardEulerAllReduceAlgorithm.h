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

    const AlgorithmType getType() override;
    const char* getTypeString() override;

    void Init() override;
    void Finalize() override;
    double Run() override;

    void StepBegin(Branch*) override;
    void StepEnd(Branch*, hpx_t) override;
    void CommStepBegin(Branch*) override;
    void CommStepEnd(Branch*) override;
    void AfterSpike(Branch*) override;

    static void SubscribeAllReduces  (hpx_t *& allReduces, size_t allReducesCount);
    static void UnsubscribeAllReduces(hpx_t *& allReduces, size_t allReducesCount);
    static void waitForSpikesDelivery(Branch * b, hpx_t spikesLco);

    const size_t allReducesCount = 1;
    static hpx_t* allReduces;

  private:
};

}; //algorithm

}; //neurox
