#pragma once
#include "neurox.h"

using namespace neurox;

namespace neurox {

namespace algorithms {

class CoreneuronAlgorithm : public Algorithm {
 public:
  CoreneuronAlgorithm();
  ~CoreneuronAlgorithm();

  const Algorithms GetId() override;
  const char* GetString() override;

  void Init() override;
  void Clear() override;
  double Launch() override;

  void StepBegin(Branch*) override;
  void StepEnd(Branch*, hpx_t) override;
  void Run(Branch*, const void*) override;
  hpx_t SendSpikes(Neuron*, double, double) override;

  class CommunicationBarrier : public AlgorithmMetadata {
   public:
    CommunicationBarrier();
    ~CommunicationBarrier();
  };
};

};  // algorithm

};  // neurox
