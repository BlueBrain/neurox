/*
# =============================================================================
# Copyright (c) 2015 - 2021-2022 Blue Brain Project/EPFL
#
# See top-level LICENSE file for details.
# =============================================================================.
*/
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

  void NeuronSyncInit(Branch*) override;
  void NeuronSyncEnd(Branch*) override;
  double GetNeuronMaxStep(Branch* b) override;
  void SendSpikes(Neuron*, double, double) override;

  class CommunicationBarrier : public SynchronizerNeuronInfo {
   public:
    CommunicationBarrier();
    ~CommunicationBarrier();
  };
};

};  // namespace synchronizers

};  // namespace neurox
