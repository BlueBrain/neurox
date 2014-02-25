#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::algorithms;

Algorithm* Algorithm::New(AlgorithmType type)
{
    switch (type)
    {
    case AlgorithmType::BackwardEulerDebugMode:
        return new BackwardEulerDebugModeAlgorithm();
    case AlgorithmType::BackwardEulerAllReduce:
        return new BackwardEulerAllReduceAlgorithm();
    case AlgorithmType::BackwardEulerSlidingTimeWindow:
        return new BackwardEulerSlidingTimeWindowAlgorithm();
    case AlgorithmType::BackwardEulerTimeDependencyLCO:
        return new BackwardEulerTimeDependencyLCOAlgorithm();
    default:
        return nullptr;
    }
    return nullptr;
}


void Algorithm::Run(AlgorithmType algorithm)
{
    int totalSteps = (inputParams->tstop - inputParams->tstart) / inputParams->dt;
    printf("neurox::Algorithm::%s (%d neurons, t=%.03f secs, dt=%.03f milisecs, %d steps)\n",
            algorithm == AlgorithmType::BackwardEulerDebugMode  ? "BackwardEulerDebugWithCommBarrier" :
           (algorithm == AlgorithmType::All                     ? "ALL" :
           (algorithm == AlgorithmType::BackwardEulerAllReduce  ? "BackwardEulerWithAllReduceBarrier" :
           (algorithm == AlgorithmType::BackwardEulerSlidingTimeWindow ? "BackwardEulerWithSlidingTimeWindow" :
           (algorithm == AlgorithmType::BackwardEulerTimeDependencyLCO ? "BackwardEulerWithTimeDependencyLCO" :
            "unknown" )))), neurons->size(), inputParams->tstop/1000, inputParams->dt, totalSteps);
    fflush(stdout);

    AlgorithmType previousAlgorithm = inputParams->algorithm;
    hpx_bcast_rsync(neurox::SetAlgorithmVariables, &algorithm, sizeof(AlgorithmType));

#ifdef NDEBUG //benchmark info
    hpx_time_t now = hpx_time_now();
#endif

    if (algorithm == AlgorithmType::BackwardEulerDebugMode)
    {
        int commStepSize = Neuron::CommunicationBarrier::commStepSize;
        for (int s=0; s<totalSteps; s+=Neuron::CommunicationBarrier::commStepSize)
        {
            #ifdef NEUROX_TIME_STEPPING_VERBOSE
              if (hpx_get_my_rank()==0)
                DebugMessage(std::string("-- t="+std::to_string(inputParams->dt*s)+" ms\n").c_str());
            #endif

            //Reduction at locality is not implemented (this mode is for debugging only)
            neurox_hpx_call_neurons_lco(Branch::BackwardEuler, &commStepSize, sizeof(int));

            #ifndef NDEBUG
              if (inputParams->parallelDataLoading) //if parallel execution... spike exchange
                hpx_bcast_rsync(neurox::input::Debugger::NrnSpikeExchange);
            #endif
        }
    }
    else
    {
        if (inputParams->allReduceAtLocality &&
            (  algorithm == AlgorithmType::BackwardEulerAllReduce
            || algorithm == AlgorithmType::BackwardEulerSlidingTimeWindow))
            hpx_bcast_rsync(Branch::BackwardEulerOnLocality, &totalSteps, sizeof(int));
        else
            neurox_hpx_call_neurons_lco(Branch::BackwardEuler, &totalSteps, sizeof(int));
    }

#ifdef NDEBUG
    //output benchmark info
    double timeElapsed = hpx_time_elapsed_ms(now)/1e3;
    printf("csv,%d,%d,%d,%.1f,%.1f,%d,%d,%d,%d,%.2f\n", neurons->size(), hpx_get_num_ranks(),
        hpx_get_num_threads(), neurons->size() / (double) hpx_get_num_ranks(), inputParams->tstop,
        algorithm, inputParams->multiMex ? 1:0, inputParams->branchingDepth,
        inputParams->allReduceAtLocality ? 1:0, timeElapsed);
    fflush(stdout);
#else
    //compare final results
    if (inputParams->branchingDepth==0)
    if (algorithm != AlgorithmType::BackwardEulerDebugMode //not fixed comm barrier
    && inputParams->parallelDataLoading) //and not serial
    {
        //re-run whole simulation and comparae final result
        DebugMessage("neurox::re-running simulation in Coreneuron to compare final result...\n");
        fflush(stdout);
        int commStepSize = Neuron::CommunicationBarrier::commStepSize;
        for (int s=0; s<totalSteps; s+=Neuron::CommunicationBarrier::commStepSize)
        {
            hpx_bcast_rsync(neurox::input::Debugger::FixedStepMinimal, &commStepSize, sizeof(int));
            hpx_bcast_rsync(neurox::input::Debugger::NrnSpikeExchange);
        }
    }
    neurox::input::Debugger::CompareAllBranches();
#endif

    hpx_bcast_rsync(neurox::SetAlgorithmVariables, &previousAlgorithm, sizeof(AlgorithmType));
}
