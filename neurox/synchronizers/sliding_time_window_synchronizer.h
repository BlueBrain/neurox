#pragma once
#include "neurox.h"

using namespace neurox;

namespace neurox {

namespace synchronizers {

class SlidingTimeWindowSynchronizer : public Synchronizer {
 public:
  SlidingTimeWindowSynchronizer();
  ~SlidingTimeWindowSynchronizer();

  const SynchronizerIds GetId() override;
  const char* GetString() override;
  void InitLocality() override;
  void ClearLocality() override;
  void BeforeSteps(Branch*) override;
  double GetMaxStep(Branch*) override;
  void AfterSteps(Branch*, hpx_t) override;
  hpx_t SendSpikes(Neuron*, double, double) override;
  double GetLocalityReductionInterval() override;
  void LocalityReduce() override;

  static const size_t kAllReducesCount = 2;

 private:
};

};  // synchronizer

};  // neurox
