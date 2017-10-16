#pragma once

#include "neurox.h"

namespace neurox {

namespace synchronizers {

/**
 * @brief The Synchronizers enum
 * Uniquely identifies the Id of each synchronizer
 */
enum class SynchronizerIds : int {
  kDebug = 0,          // For debug only
  kTimeDependency = 1, /* Needs to be first */
  kAllReduce = 2,
  kSlidingTimeWindow = 3,
  kCoreneuron = 4,
  kSynchronizersCount = 5,
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

  /// Initialize synchronizer meta data in locality
  virtual void InitLocality() {}

  /// Initialize synchronizer meta data in neuron
  virtual void InitNeuron(Branch*) {}

  /// Clears/finalizes synchronizer meta data at locality-level
  virtual void ClearLocality() {}

  /// Clears/finalizes synchronizer meta data at neuron-level
  virtual void ClearNeuron(Branch*) {}

  /// To be called at beginning of steps
  virtual void BeforeSteps(Branch*) {}

  /// Maximum possible step, without synchronization
  virtual double GetMaxStep(Branch*) = 0;

  /// To be called at end of steps
  virtual void AfterSteps(Branch*, hpx_t spikesLco) {}

  /// To handle sending of spikes
  virtual hpx_t SendSpikes(Neuron* n, double tt, double t) {
    return Neuron::SendSpikesAsync(n, tt);
  }

  /// To handle receival of spikes
  virtual void AfterReceiveSpikes(Branch*, hpx_t, neuron_id_t, spike_time_t,
                                  spike_time_t) {}

  /// Time-step between locality-based comm. reductions
  /// (default 0, filtered as 'no locality reduction')
  virtual double GetLocalityReductionInterval() { return 0.; }

  /// Locatility-based reduction, at every reduction-interval
  virtual void LocalityReduce() {}

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
