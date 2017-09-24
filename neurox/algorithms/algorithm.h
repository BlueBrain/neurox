#pragma once

#include "neurox.h"

namespace neurox {

namespace algorithms {

/**
 * @brief The AlgorithmIds enum
 * Uniquely identifies the Id of each algorithm
 */
enum class AlgorithmId : int {
  kBackwardEulerDebug = 0,  // For debug only
  kBackwardEulerAllReduce = 1,
  kBackwardEulerSlidingTimeWindow = 2,
  kBackwardEulerTimeDependencyLCO = 3,
  kBackwardEulerCoreneuron = 4,
  kCvodes = 5,
  kBenchmarkAll = 9  // Benchmark of all non-debug modes
};

/**
 * @brief The NeuronMetaData class
 * Purely abstract class, represents metadata at neuron level
 * relative to a neuron
 */
class AlgorithmMetadata{
  public:
    /// Returns an instantiated metadata for the algorithm of given type
    static AlgorithmMetadata* New(AlgorithmId);
};

/**
 * @brief The Algorithm class
 * Virtual class describing all methods required to be implemented by all
 * different algorithms
 */
class Algorithm {
 public:

  Algorithm(){};
  virtual ~Algorithm(){};

  /// Init for Fixed Step methods
  static void FixedStepMethodsInit();

  /// Returns total count of steps of fixed size defined by user
  static int GetTotalStepsCount();

  /// Returns an instantiated class of the given type
  static Algorithm* New(AlgorithmId);

  /// Prints info on algorithm, data size, and steps count
  void PrintStartInfo();

  /// Returns class type
  const virtual AlgorithmId GetId() = 0;

  /// Returns class type as string
  const virtual char* GetString() = 0;

  /// Launch simulation on all neurons or localities
  virtual double Launch() = 0;

  /// Runs simulation for given branch
  virtual void Run(Branch*, const void*) = 0;

  /// Initialize algorithm meta data
  virtual void Init(){};

  /// Clears/finalizes algorithm meta data
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

};  // algorithms

};  // neurox

// TODO can we move this somewhere else?
#include "neurox/algorithms/allreduce_algorithm.h"
#include "neurox/algorithms/coreneuron_algorithm.h"
#include "neurox/algorithms/debug_algorithm.h"
#include "neurox/algorithms/sliding_time_window_algorithm.h"
#include "neurox/algorithms/time_dependency_lco_algorithm.h"
#include "neurox/algorithms/cvodes_algorithm.h"
