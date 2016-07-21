#pragma once
#include "neurox.h"

using namespace neurox;

namespace neurox {

namespace algorithms {

class AllReduceAlgorithm : public Algorithm {
 public:
  AllReduceAlgorithm();
  ~AllReduceAlgorithm();

  const AlgorithmType GetType() override;
  const char* GetTypeString() override;

  void Init() override;
  void Clear() override;
  double Launch() override;

  void StepBegin(Branch*) override;
  void StepEnd(Branch*, hpx_t) override;
  void Run(Branch*, const void*) override;
  hpx_t SendSpikes(Neuron*, double, double) override;

  static void SubscribeAllReduces(hpx_t*& allReduces, size_t allReducesCount);
  static void UnsubscribeAllReduces(hpx_t*& allReduces, size_t allReducesCount);
  static void WaitForSpikesDelivery(Branch* b, hpx_t spikesLco);
  static hpx_t SendSpikes2(Neuron*, double);
  static void Run2(Branch*, const void*);

  const size_t allReducesCount = 1;
  static hpx_t* allReduces;

  class AllReducesInfo : public AlgorithmMetaData {
   public:
    AllReducesInfo();
    ~AllReducesInfo();

    static void RegisterHpxActions();  ///> Register all HPX actions

    // set by initNodeLevelInformaion
    static int reductionsPerCommStep;

    // for node level reduction only (initialized by initNodeLevelInformaion)
    class AllReduceLocality {
     public:
      static std::vector<hpx_t>* localityNeurons;
      static hpx_t* allReduceFuture;
      static hpx_t* allReduceLco;
      static int* allReduceId;

      static hpx_action_t SubscribeAllReduce;
      static hpx_action_t UnsubscribeAllReduce;

      static int SubscribeAllReduce_handler(const hpx_t*, const size_t);
      static int UnsubscribeAllReduce_handler(const hpx_t*, const size_t);
    };

    // initiated by constructor (one per neuron)
    std::queue<hpx_t> spikesLcoQueue;
    hpx_t* allReduceFuture;
    hpx_t* allReduceLco;
    int* allReduceId;

    // all-reduce functions
    static hpx_action_t Init;
    static hpx_action_t Reduce;
    static hpx_action_t SubscribeAllReduce;
    static hpx_action_t UnsubscribeAllReduce;
    static hpx_action_t SetReductionsPerCommStep;

    static void Init_handler(const void*, const size_t);
    static void Reduce_handler(void* rhs, const void* lhs, const size_t);
    static int SubscribeAllReduce_handler(const hpx_t*, const size_t);
    static int UnsubscribeAllReduce_handler(const hpx_t*, const size_t);
    static int SetReductionsPerCommStep_handler(const int*, const size_t);
  };

 private:
};

};  // algorithm

};  // neurox
