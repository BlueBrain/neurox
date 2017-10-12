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

  /// Init for Fixed Step methods
  static void FixedStepMethodsInit();

  /// Returns an instantiated class of the given type
  static Synchronizer* New(Synchronizers);

  /// Prints info on synchronizer, data size, and steps count
  void PrintStartInfo();

  /// Returns class type
  const virtual Synchronizers GetId() = 0;

  /// Returns class type as string
  const virtual char* GetString() = 0;

  /// Launch simulation on all neurons or localities
  virtual void Launch() = 0;

  /// Runs simulation for given branch
  virtual void Run(Branch*, const void*) = 0;

  /// Initialize synchronizer meta data
  virtual void Init(){};

  /// Clears/finalizes synchronizer meta data
  virtual void Clear(){};

  /// To be called at beginning of step
  virtual void StepBegin(Branch*){};

  /// To be called at end of step
  virtual void StepEnd(Branch*, hpx_t spikesLco){};

  /// To handle sending of spikes
  virtual hpx_t SendSpikes(Neuron*, double tt, double t) = 0;

  /// To handle receival of spikes
  virtual void AfterReceiveSpikes(Branch*, hpx_t, neuron_id_t, spike_time_t,
                                  spike_time_t){};

  static hpx_action_t Init;
  static hpx_action_t Clear;
  static hpx_action_t RunNeuron;
  static hpx_action_t RunLocality;

  static void RegisterHpxActions();  ///> Register all HPX actions

private:

  static int Init_handler(const int*, const size_t);
  static int Clear_handler();
  static int RunNeuron_handler(const double*, const size_t);
  static int RunLocality_handler(const double*, const size_t);

};

};  // synchronizers

};  // neurox
