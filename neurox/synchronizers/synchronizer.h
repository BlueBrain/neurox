#pragma once

#include "neurox.h"

namespace neurox {

namespace synchronizers {

/**
 * @brief The Synchronizers enum
 * Uniquely identifies the Id of each synchronizer
 */
enum class SynchronizerIds : int {
  kTimeDependency = 0, /* Needs to be first */
  kAllReduce = 1,
  kSlidingTimeWindow = 2,
  kCoreneuron = 3,
  kSynchronizersCount = 3,
  kDebug = 8,        // For debug only
  kBenchmarkAll = 9  // Benchmark of all non-debug modes
};

/**
 * @brief The SynchronizerMetadata class
 * Purely abstract class, represents metadata at neuron level
 * relative to a neuron
 */
class SynchronizerNeuronInfo {
 public:
  /// Returns an instantiated metadata for the synchronizer of given type
  static SynchronizerNeuronInfo* New(SynchronizerIds);
};

/**
 * @brief The SyncSynchronizer class
 * Virtual class describing all methods required to be implemented by all
 * different synchronizers
 */
class Synchronizer {
 public:
  Synchronizer(){};
  virtual ~Synchronizer(){};

  /// Returns an instantiated class of the given type
  static Synchronizer* New(SynchronizerIds);

  /// Returns class type
  const virtual SynchronizerIds GetId() = 0;

  /// Returns class type as string
  const virtual char* GetString() = 0;




  /// Initialize synchronizer meta data in neuron
  virtual void InitNeuron(Branch*) {}

  /// Clears/finalizes synchronizer meta data at neuron-level
  virtual void ClearNeuron(Branch*) {}

  /// neuron-based reduction: interval between reduction
  virtual double NeuronSyncInterval(Branch*) = 0;

  /// Neuron-based reduction: called at the start of every neuron reduction
  virtual void NeuronSyncInit(Branch* b) {assert(b->soma_);}

  /// Neuron-based reduction: called at the end of every neuron reduction
  virtual void NeuronSyncEnd(Branch* b, hpx_t spikesLco) {assert(b->soma_);}




  /// Initialize synchronizer meta data in locality
  virtual void InitLocality() {}

  /// Clears/finalizes synchronizer meta data at locality-level
  virtual void ClearLocality() {}

  /* Locatility-based reduction:
   * Time-step between locality-based comm. reductions:
   *  - positive value: call LocalityReduce at every interval
   *  - 0 (default): no reduction, launch neurons independently
   *  - -1 : locality level syncrhonizer, launch last neuron first
   */
  virtual double LocalitySyncInterval() { return 0; }

  /// Locatility-based reduction: called at the start every reduction-interval
  virtual void LocalitySyncInit() {}

  /// Locatility-based reduction: called at the end of every reduction-interval
  virtual void LocalitySyncEnd() {}



  /// spikes handling: how it hadles outgoing spikes after Action-Potential
  virtual hpx_t SendSpikes(Neuron* n, double tt, double t) = 0;

  /// spikes handling: how it reacts to receival of spikes
  virtual void AfterReceiveSpikes(Branch*, hpx_t, neuron_id_t, spike_time_t,
                                  spike_time_t) {}


  static hpx_action_t CallInitLocality;
  static hpx_action_t CallInitNeuron;
  static hpx_action_t RunNeuron;
  static hpx_action_t RunLocality;
  static hpx_action_t CallClearLocality;
  static hpx_action_t CallClearNeuron;
  static hpx_action_t NeuronInfoConstructor;
  static hpx_action_t NeuronInfoDestructor;

  static void RegisterHpxActions();  ///> Register all HPX actions

  /// auxiliar method for CallLocalNeurons
  /// (TODO move it to wrappers?)
  static hpx_action_t CallAllNeuronsAux;

 private:
  static int CallInitLocality_handler(const int*, const size_t);
  static int CallInitNeuron_handler();
  static int RunNeuron_handler(const double*, const size_t);
  static int RunLocality_handler(const double*, const size_t);
  static int CallClearLocality_handler();
  static int CallClearNeuron_handler();
  static int NeuronInfoConstructor_handler(const int*, const size_t);
  static int NeuronInfoDestructor_handler();

  static int CallAllNeuronsAux_handler(const int, const void* [],
                                       const size_t[]);
};

};  // synchronizers

};  // neurox
