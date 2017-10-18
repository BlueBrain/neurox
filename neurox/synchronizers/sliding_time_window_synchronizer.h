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
  void NeuronReduceInit(Branch*) override;
  hpx_t SendSpikes(Neuron*, double, double) override;
  void NeuronReduceEnd(Branch*, hpx_t) override;
  double GetReduceInterval() override;
  void LocalityReduceInit() override;

  static const size_t kAllReducesCount = 2;

 private:
};

};  // synchronizer

};  // neurox
