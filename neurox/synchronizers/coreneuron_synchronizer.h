#pragma once
#include "neurox.h"

using namespace neurox;

namespace neurox {

namespace synchronizers {

class CoreneuronSynchronizer : public Synchronizer {
 public:
  CoreneuronSynchronizer();
  ~CoreneuronSynchronizer();

  const SynchronizerIds GetId() override;
  const char* GetString() override;

  void InitLocality() override;
  void ClearLocality() override;

  void BeforeSteps(Branch*) override;
  void AfterSteps(Branch*, hpx_t) override;
  double GetMaxStep(Branch* b) override;
  hpx_t SendSpikes(Neuron*, double, double) override;

  class CommunicationBarrier : public SynchronizerNeuronInfo {
   public:
    CommunicationBarrier();
    ~CommunicationBarrier();
  };
};

};  // synchronizer

};  // neurox
