#pragma once

#include "neurox/neurox.h"

namespace neurox
{

namespace algorithms
{

enum AlgorithmType {
    BackwardEulerCoreneuronDebug=-1, //For debug only
    BackwardEulerCoreneuron=0,
    BackwardEulerAllReduce=1,
    BackwardEulerSlidingTimeWindow=2,
    BackwardEulerTimeDependencyLCO=3,
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
