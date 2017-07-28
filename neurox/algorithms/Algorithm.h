#pragma once
#include "neurox.h"
#include "neurox/algorithms/AlgorithmMetaData.h"

namespace neurox {

namespace algorithms {

/**
 * @brief The Algorithm class
 * Virtual class describing all methods required to be implemented by all
 * different algorithms
 */
class Algorithm {
 public:
  Algorithm(){};
  virtual ~Algorithm(){};

  /// Returns total count of steps of fixed size defined by user
  static int getTotalStepsCount();

  /// Returns an instantiated class of the given type
  static Algorithm* New(AlgorithmType);

  /// Prints info on algorithm, data size, and steps count
  void PrintStartInfo();

  /// Returns class type
  const virtual AlgorithmType getType() = 0;

  /// Returns class type as string
  const virtual char* getTypeString() = 0;

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

#include "neurox/algorithms/AllReduceAlgorithm.h"
#include "neurox/algorithms/CoreneuronAlgorithm.h"
#include "neurox/algorithms/DebugAlgorithm.h"
#include "neurox/algorithms/SlidingTimeWindowAlgorithm.h"
#include "neurox/algorithms/TimeDependencyLCOAlgorithm.h"
