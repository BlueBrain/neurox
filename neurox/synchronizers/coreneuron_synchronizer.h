#pragma once
#include "neurox.h"

using namespace neurox;

namespace neurox {

namespace synchronizers {

class CoreneuronSynchronizer : public Synchronizer {
 public:
  CoreneuronSynchronizer();
  ~CoreneuronSynchronizer();

  const Synchronizers GetId() override;
  const char* GetString() override;

  void Init() override;
  void Clear() override;

  void BeforeStep(Branch*) override;
  void AfterStep(Branch*, hpx_t) override;
  double GetMaxStepTime(Branch* b) override;
  hpx_t SendSpikes(Neuron*, double, double) override;

  class CommunicationBarrier : public SynchronizerMetadata {
   public:
    CommunicationBarrier();
    ~CommunicationBarrier();
  };
};

};  // synchronizer

};  // neurox
