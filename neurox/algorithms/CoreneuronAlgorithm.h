#pragma once
#include "neurox.h"

using namespace neurox;

namespace neurox {

namespace algorithms {

class CoreneuronAlgorithm : public Algorithm {
 public:
  CoreneuronAlgorithm();
  ~CoreneuronAlgorithm();

  const AlgorithmType getType() override;
  const char* getTypeString() override;

  void Init() override;
  void Clear() override;
  double Launch() override;

  void StepBegin(Branch*) override;
  void StepEnd(Branch*, hpx_t) override;
  void Run(Branch*, const void*) override;
  hpx_t SendSpikes(Neuron*, double, double) override;

  class CommunicationBarrier : public AlgorithmMetaData {
   public:
    CommunicationBarrier();
    ~CommunicationBarrier();

    static constexpr int commStepSize = 4;  ///> Fixed communication step size
  };
};

};  // algorithm

};  // neurox
