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
  kTimeDependencyLCO = 3,
  kCoreneuron = 4,
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

  /// Returns total count of steps of fixed size defined by user
  static int GetTotalStepsCount();

  /// Returns an instantiated class of the given type
  static Synchronizer* New(Synchronizers);

  /// Prints info on synchronizer, data size, and steps count
  void PrintStartInfo();

  /// Returns class type
  const virtual Synchronizers GetId() = 0;

  /// Returns class type as string
  const virtual char* GetString() = 0;

  /// Launch simulation on all neurons or localities
  virtual double Launch() = 0;

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
};

};  // synchronizers

};  // neurox

// TODO can we move this somewhere else?
#include "neurox/synchronizers/allreduce_synchronizer.h"
#include "neurox/synchronizers/coreneuron_synchronizer.h"
#include "neurox/synchronizers/debug_synchronizer.h"
#include "neurox/synchronizers/sliding_time_window_synchronizer.h"
#include "neurox/synchronizers/time_dependency_lco_synchronizer.h"
