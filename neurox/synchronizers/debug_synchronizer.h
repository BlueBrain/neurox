#pragma once
#include "neurox.h"

using namespace neurox;

namespace neurox {

namespace synchronizers {

class DebugSynchronizer : public Synchronizer {
 public:
  DebugSynchronizer();
  ~DebugSynchronizer();

  const SynchronizerIds GetId() override;
  const char* GetString() override;

  void InitLocality() override;
  void Clear() override;

  void BeforeStep(Branch*) override;
  void AfterStep(Branch*, hpx_t) override;
  double GetMaxStepTime(Branch*) override;
  hpx_t SendSpikes(Neuron*, double, double) override;

  void Launch();
  void Run(Branch*, const void*);

  class CommunicationBarrier : public SynchronizerNeuronInfo {
   public:
    CommunicationBarrier();
    ~CommunicationBarrier();

    /// LCO for all spikes of previous Comm Step (for fixed step methods)
    hpx_t all_spikes_lco_;
  };
};

};  // synchronizer

};  // neurox
