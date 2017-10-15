#pragma once
#include "neurox.h"

using namespace neurox;

namespace neurox {

namespace synchronizers {

class AllreduceSynchronizer : public Synchronizer {
 public:
  AllreduceSynchronizer();
  ~AllreduceSynchronizer();

  const SynchronizerIds GetId() override;
  const char* GetString() override;
  void InitLocality() override;
  void InitNeuron(Branch *b) override;
  void ClearLocality() override;
  void ClearNeuron(Branch*) override;
  void BeforeSteps(Branch*) override;
  double GetMaxStepTime(Branch*) override;
  void AfterSteps(Branch*, hpx_t) override;
  hpx_t SendSpikes(Neuron*, double, double) override;
  double GetLocalityReductionInterval() override;
  void LocalityReduce() override;

  static void SubscribeAllReducesLocality(size_t allreduces_count);
  static void SubscribeAllReducesNeuron(Branch *b, size_t allreduces_count);
  static void UnsubscribeAllReducesLocality(size_t allreduces_count);
  static void UnsubscribeAllReducesNeuron(Branch* b, size_t allreduces_count);
  static void WaitForSpikesDelivery(Branch* b, hpx_t spikes_lco);

  static const size_t kAllReducesCount = 1;
  static hpx_t* allreduces_;

  static void NeuronReduce(const Branch*, const int);
  static double GetMaxStepTime2(const Branch*, const int);
  static double GetLocalityReductionInterval2(const double);

  // for node level reduction only (initialized by initNodeLevelInformaion)
  class AllReduceLocalityInfo {
   public:
    static hpx_t* allreduce_future_;
    static hpx_t* allreduce_lco_;
    static int* allreduce_id_;
    static int next_allreduce_id;

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

    static void RegisterHpxActions();  ///> Register all HPX actions

    // initiated by constructor (one per neuron)
    std::queue<hpx_t> spikes_lco_queue_;
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

};  // synchronizer

};  // neurox
