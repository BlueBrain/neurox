#pragma once

#include "neurox/neurox.h"

namespace neurox {

namespace algorithms {

enum AlgorithmType {
  BackwardEulerDebug = 0,  // For debug only
  BackwardEulerAllReduce = 1,
  BackwardEulerSlidingTimeWindow = 2,
  BackwardEulerTimeDependencyLCO = 3,
  BackwardEulerCoreneuron = 4,
  BenchmarkEnd = 3,
  All = 9  // Benchmark of all non-debug modes
};

class AlgorithmMetaData {
 public:
  static AlgorithmMetaData* New(AlgorithmType);
};

};  // algorithms

};  // neurox
