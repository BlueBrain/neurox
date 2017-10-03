#pragma once
#include "neurox.h"

using namespace neurox;

namespace neurox {

namespace algorithms {

class SlidingTimeWindowAlgorithm : public Algorithm {
 public:
  SlidingTimeWindowAlgorithm();
  ~SlidingTimeWindowAlgorithm();

  const SyncAlgorithms GetId() override;
  const char* GetString() override;

  void Init() override;
  void Clear() override;
  double Launch() override;

  void StepBegin(Branch*) override;
  void StepEnd(Branch*, hpx_t) override;
  void Run(Branch*, const void*) override;
  hpx_t SendSpikes(Neuron*, double, double) override;

  const size_t kAllReducesCount = 2;
  static hpx_t* allreduces_;

 private:
};

};  // algorithm

};  // neurox
