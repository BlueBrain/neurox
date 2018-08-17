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
  void NeuronSyncInit(Branch*) override;
  hpx_t SendSpikes(Neuron*, double, double) override;
  double GetNeuronMaxStep(Branch*) override;
  void NeuronSyncEnd(Branch*, hpx_t) override;
  double LocalitySyncInterval() override;
  void LocalitySyncInit() override;

  static const size_t kAllReducesCount = 2;

 private:
};

};  // namespace synchronizers

};  // namespace neurox
