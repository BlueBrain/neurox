#pragma once

#include "neurox/neurox.h"

namespace neurox
{

namespace algorithms
{

enum AlgorithmType {
    BackwardEulerCoreneuronDebug=-1, //For debug only
    BackwardEulerAllReduce=0,
    BackwardEulerSlidingTimeWindow=1,
    BackwardEulerTimeDependencyLCO=2,
    BenchmarkEnd=3,
    All=9 //Benchmark of all non-debug modes
};

class AlgorithmMetaData
{
  public:
    static AlgorithmMetaData* New(AlgorithmType);
};

}; //algorithms

}; //neurox
