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
  double Launch() override;

  void StepBegin(Branch*) override;
  void StepEnd(Branch*, hpx_t) override;
  void Run(Branch*, const void*) override;
  hpx_t SendSpikes(Neuron*, double, double) override;

  static void SubscribeAllReduces(hpx_t*& allreduces, size_t allreduces_count);
  static void UnsubscribeAllReduces(hpx_t*& allreduces,
                                    size_t allreduces_count);
  static void WaitForSpikesDelivery(Branch* b, hpx_t spikes_lco);
  static void Run2(Branch*, const void*);

  const size_t kAllReducesCount = 1;
  static hpx_t* allreduces_;

  class AllReducesInfo : public SynchronizerMetadata {
   public:
    AllReducesInfo();
    ~AllReducesInfo();

    static void RegisterHpxActions();  ///> Register all HPX actions

    // set by initNodeLevelInformaion
    static int reductions_per_comm_step_;

    // for node level reduction only (initialized by initNodeLevelInformaion)
    class AllReduceLocality {
     public:
      static std::vector<hpx_t>* locality_neurons_;
      static hpx_t* allreduce_future_;
      static hpx_t* allreduce_lco_;
      static int* allreduce_id_;

      static hpx_action_t SubscribeAllReduce;
      static hpx_action_t UnsubscribeAllReduce;

      static int SubscribeAllReduce_handler(const hpx_t*, const size_t);
      static int UnsubscribeAllReduce_handler(const hpx_t*, const size_t);
    };

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
    static hpx_action_t SetReductionsPerCommStep;

    static void Init_handler(void*, const size_t);
    static void Reduce_handler(void* rhs, const void* lhs, const size_t);
    static int SubscribeAllReduce_handler(const hpx_t*, const size_t);
    static int UnsubscribeAllReduce_handler(const hpx_t*, const size_t);
    static int SetReductionsPerCommStep_handler(const int*, const size_t);
  };

 private:
};

};  // synchronizer

};  // neurox
