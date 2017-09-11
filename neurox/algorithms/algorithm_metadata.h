#pragma once

#include "neurox/neurox.h"

namespace neurox {

namespace algorithms {

enum AlgorithmType {
  kBackwardEulerDebug = 0,  // For debug only
  kBackwardEulerAllReduce = 1,
  kBackwardEulerSlidingTimeWindow = 2,
  kBackwardEulerTimeDependencyLCO = 3,
  kBackwardEulerCoreneuron = 4,
  kCvodes = 5,
  kBenchmarkAll = 9  // Benchmark of all non-debug modes
};

class AlgorithmMetadata {
 public:
  static AlgorithmMetadata* New(AlgorithmType);
};

};  // algorithms

};  // neurox
