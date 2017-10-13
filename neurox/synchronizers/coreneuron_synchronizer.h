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

  void StepBegin(Branch*) override;
  void StepEnd(Branch*, hpx_t) override;
  void Run(Branch*, const void*) override;
  hpx_t SendSpikes(Neuron*, double, double) override;

  void Launch();

  class CommunicationBarrier : public SynchronizerMetadata {
   public:
    CommunicationBarrier();
    ~CommunicationBarrier();
  };
};

};  // synchronizer

};  // neurox
