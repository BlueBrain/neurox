#pragma once
#include "neurox.h"

using namespace neurox;

namespace neurox {

namespace synchronizers {

class AllreduceSynchronizer : public Synchronizer {
 public:
  AllreduceSynchronizer();
  ~AllreduceSynchronizer();

  const Synchronizers GetId() override;
  const char* GetString() override;

  void Init() override;
  void Clear() override;
  void BeforeStep(Branch*) override;
  void AfterStep(Branch*, hpx_t) override;
  hpx_t SendSpikes(Neuron*, double, double) override;
  double GetLocalityReductionInterval() override;
  void LocalityReduce() override;

  void Run(Branch*, const void*);

  static void SubscribeAllReduces(size_t allreduces_count);
  static void UnsubscribeAllReduces(size_t allreduces_count);
  static void WaitForSpikesDelivery(Branch* b, hpx_t spikes_lco);
  static void Run2(Branch*, const void*);

  const size_t kAllReducesCount = 1;
  static hpx_t* allreduces_;

  // for node level reduction only (initialized by initNodeLevelInformaion)
  class AllReduceLocalityInfo {
   public:
    static hpx_t* allreduce_future_;
    static hpx_t* allreduce_lco_;
    static int* allreduce_id_;
    static int next_allreduce_id;

    static void LocalityReduce(int);

    static hpx_action_t SubscribeAllReduce;
    static hpx_action_t UnsubscribeAllReduce;

    static int SubscribeAllReduce_handler(const hpx_t*, const size_t);
    static int UnsubscribeAllReduce_handler(const hpx_t*, const size_t);
  };

  class AllReducesInfo : public SynchronizerMetadata {
   public:
    AllReducesInfo();
    ~AllReducesInfo();

    static void RegisterHpxActions();  ///> Register all HPX actions

    // initiated by constructor (one per neuron)
    std::queue<hpx_t> spikes_lco_queue_;
    hpx_t* allreduce_future_;
    hpx_t* allreduce_lco_;
    int* allreduce_id_;

    // all-reduce functions
    static hpx_action_t Init;
    static hpx_action_t Reduce;
    static hpx_action_t SubscribeAllReduce;
    static hpx_action_t UnsubscribeAllReduce;

    static void Init_handler(void*, const size_t);
    static void Reduce_handler(void* rhs, const void* lhs, const size_t);
    static int SubscribeAllReduce_handler(const hpx_t*, const size_t);
    static int UnsubscribeAllReduce_handler(const hpx_t*, const size_t);
  };

 private:
};

};  // synchronizer

};  // neurox
