/*
# =============================================================================
# Copyright (c) 2015 - 2021 Blue Brain Project/EPFL
#
# See top-level LICENSE file for details.
# =============================================================================.
*/
#pragma once
#include "neurox.h"

using namespace neurox;

namespace neurox {

namespace synchronizers {

typedef std::pair<floble_t, hpx_t> TimedSpike;

class AllreduceSynchronizer : public Synchronizer {
 public:
  AllreduceSynchronizer();
  ~AllreduceSynchronizer();

  const SynchronizerIds GetId() override;
  const char* GetString() override;
  void InitLocality() override;
  void ClearLocality() override;
  void NeuronSyncInit(Branch*) override;
  void SendSpikes(Neuron*, double, double) override;
  double GetNeuronMaxStep(Branch*) override;
  void NeuronSyncEnd(Branch*) override;
  double LocalitySyncInterval() override;
  void LocalitySyncInit() override;

  static void SubscribeAllReduces(size_t allreduces_count);
  static void UnsubscribeAllReduces(size_t allreduces_count);
  static void WaitForSpikesDelivery(Branch* b);

  static const size_t kAllReducesCount = 1;
  static hpx_t* allreduces_;

  static void NeuronReduce(const Branch*, const int);
  static double NeuronReduceInterval2(const int);
  static double LocalityReduceInterval2(const double);
  static void SendSpikes2(Neuron*, spike_time_t t);

  static void RegisterHpxActions();  ///> Register all HPX actions

  // for node level reduction only (initialized by initNodeLevelInformaion)
  class AllReduceLocalityInfo {
   public:
    static hpx_t* allreduce_future_;
    static hpx_t* allreduce_lco_;
    static int* allreduce_id_;
    static int next_allreduce_id_;

    static void LocalityReduce(int);

    static hpx_action_t Subscribe;
    static hpx_action_t Unsubscribe;

    static int Subscribe_handler(const hpx_t*, const size_t);
    static int Unsubscribe_handler(const hpx_t*, const size_t);
  };

  class AllReduceNeuronInfo : public SynchronizerNeuronInfo {
   public:
    AllReduceNeuronInfo() = delete;
    AllReduceNeuronInfo(const size_t);
    ~AllReduceNeuronInfo();

    // initiated by constructor (one per neuron)
    // (At time X wait for HPX synapses Y)
    std::priority_queue<TimedSpike, std::vector<TimedSpike>,
                        std::greater_equal<TimedSpike> >
        spikes_lco_queue_;
    hpx_t* allreduce_future_;
    hpx_t* allreduce_lco_;
    int* allreduce_id_;
    int next_allreduce_id_;

    // all-reduce functions
    static hpx_action_t Init;
    static hpx_action_t Reduce;
    static hpx_action_t Subscribe;
    static hpx_action_t Unsubscribe;

    static void Init_handler(void*, const size_t);
    static void Reduce_handler(void* rhs, const void* lhs, const size_t);
    static int Subscribe_handler(const hpx_t*, const size_t);
    static int Unsubscribe_handler(const hpx_t*, const size_t);
  };

 private:
};

};  // namespace synchronizers

};  // namespace neurox
