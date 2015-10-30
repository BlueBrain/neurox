#pragma once
#include "neurox.h"

namespace neurox
{

namespace algorithms
{

enum AlgorithmType {
    BackwardEulerDebugMode=-1, //For debug only
    BackwardEulerAllReduce=0,
    BackwardEulerSlidingTimeWindow=1,
    BackwardEulerTimeDependencyLCO=2,
    BenchmarkEnd=3,
    All=9 //Benchmark of all non-debug modes
};

/**
 * @brief The Algorithm class
 * Virtual class describing all methods required to be implemented by all different algorithms
 */
class Algorithm
{
  public:
    Algorithm() {};
    virtual ~Algorithm() {};

    static int getTotalStepsCount();

    /// Returns an instantiated class of the given type
    static Algorithm* New(AlgorithmType);

    /// Prints info on algorithm, data size, and steps count
    void PrintStartInfo();

    const virtual AlgorithmType getType() = 0; ///> Returns class type
    const virtual char* getTypeString() = 0; ///> Returns class type string

    virtual void Init() {};
    virtual void Clear() {};
    virtual double Launch() = 0;
    virtual void SimulationBegin(Branch *) {};
    virtual void SimulationEnd(Branch *) {};
    virtual void StepBegin(Branch*) {};
    virtual void StepEnd(Branch*, hpx_t spikesLco) {};
    virtual void CommStepBegin(Branch*) {};
    virtual void CommStepEnd(Branch*) {};
    virtual hpx_t SendSpikes(Neuron*, double tt, double t) = 0;
    virtual void afterSpikeReceival(Branch *, hpx_t, neuron_id_t, spike_time_t, spike_time_t) {};
};

}; //algorithms

}; //neurox

#include "neurox/algorithms/BackwardEulerAllReduceAlgorithm.h"
#include "neurox/algorithms/BackwardEulerDebugModeAlgorithm.h"
#include "neurox/algorithms/BackwardEulerSlidingTimeWindowAlgorithm.h"
#include "neurox/algorithms/BackwardEulerTimeDependencyLCOAlgorithm.h"
