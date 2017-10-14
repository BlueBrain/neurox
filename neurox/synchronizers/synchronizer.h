#pragma once

#include "neurox.h"

namespace neurox {

namespace synchronizers {

/**
 * @brief The Synchronizers enum
 * Uniquely identifies the Id of each synchronizer
 */
enum class Synchronizers : int {
  kDebug = 0,  // For debug only
  kAllReduce = 1,
  kSlidingTimeWindow = 2,
  kTimeDependency = 3,
  kCoreneuron = 4,
  kSynchronizersCount = 5,
  kBenchmarkAll = 9  // Benchmark of all non-debug modes
};

/**
 * @brief The SynchronizerMetadata class
 * Purely abstract class, represents metadata at neuron level
 * relative to a neuron
 */
class SynchronizerMetadata {
 public:
  /// Returns an instantiated metadata for the synchronizer of given type
  static SynchronizerMetadata* New(Synchronizers);
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
  static Synchronizer* New(Synchronizers);

  /// Returns class type
  const virtual Synchronizers GetId() = 0;

  /// Returns class type as string
  const virtual char* GetString() = 0;

  /// Initialize synchronizer meta data in locality
  virtual void InitLocality() {}

  /// Initialize synchronizer meta data in neuron
  virtual void InitNeuron(Branch*) {}

  /// Clears/finalizes synchronizer meta data
  virtual void Clear() {}

  /// To be called at beginning of step
  virtual void BeforeStep(Branch*) {}

  virtual double GetMaxStepTime(Branch*) = 0;

  /// To be called at end of step
  virtual void AfterStep(Branch*, hpx_t spikesLco) {}

  /// To handle sending of spikes
  virtual hpx_t SendSpikes(Neuron*, double tt, double t) = 0;

  /// To handle receival of spikes
  virtual void AfterReceiveSpikes(Branch*, hpx_t, neuron_id_t, spike_time_t,
                                  spike_time_t) {}

  /// Time-step between locality-based comm. reductions
  virtual double GetLocalityReductionInterval();

  /// Locatility-based reduction
  virtual void LocalityReduce() {}

  static hpx_action_t InitializeLocality;
  static hpx_action_t InitializeNeuron;
  static hpx_action_t ClearLocality;
  static hpx_action_t RunNeuron;
  static hpx_action_t RunLocality;

  static void RegisterHpxActions();  ///> Register all HPX actions

 private:
  ///  hpx address of all neurons in this locality
  static hpx_t* locality_neurons_;

  /// length of locality_neuronx_
  static int locality_neurons_count_;

  static int InitializeLocality_handler(const int*, const size_t);
  static int InitializeNeuron_handler();
  static int ClearLocality_handler();
  static int RunNeuron_handler(const double*, const size_t);
  static int RunLocality_handler(const double*, const size_t);
};

};  // synchronizers

};  // neurox
