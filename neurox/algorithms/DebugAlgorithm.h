#pragma once
#include "neurox.h"

using namespace neurox;

namespace neurox
{

namespace algorithms
{

class DebugAlgorithm : public Algorithm
{
  public:
    DebugAlgorithm();
    ~DebugAlgorithm();

    const AlgorithmType getType() override;
    const char* getTypeString() override;

    void Init() override;
    void Clear() override;
    double Launch() override;

    void StepBegin(Branch*) override;
    void StepEnd(Branch*, hpx_t) override;
    void Run(Branch*, const void*) override;
    hpx_t SendSpikes(Neuron*, double, double) override;

    class CommunicationBarrier : public AlgorithmMetaData
    {
      public:
        CommunicationBarrier();
        ~CommunicationBarrier();

        hpx_t allSpikesLco; ///> LCO for all spikes of previous Comm Step (for fixed step methods and debug)
    };
};

}; //algorithm

}; //neurox
