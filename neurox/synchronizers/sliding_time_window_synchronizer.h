#pragma once
#include "neurox.h"

using namespace neurox;

namespace neurox {

namespace synchronizers {

class SlidingTimeWindowSynchronizer : public Synchronizer {
 public:
  SlidingTimeWindowSynchronizer();
  ~SlidingTimeWindowSynchronizer();

  const Synchronizers GetId() override;
  const char* GetString() override;

  void Init() override;
  void Clear() override;

  void BeforeStep(Branch*) override;
  void AfterStep(Branch*, hpx_t) override;
  void Run(Branch*, const void*) override;
  hpx_t SendSpikes(Neuron*, double, double) override;

  void Launch();

  const size_t kAllReducesCount = 2;
  static hpx_t* allreduces_;

 private:
};

};  // synchronizer

};  // neurox
