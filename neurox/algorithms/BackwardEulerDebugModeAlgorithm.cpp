#include "neurox/algorithms/BackwardEulerDebugModeAlgorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

DERIVED_CLASS_NAME::DERIVED_CLASS_NAME() {}

DERIVED_CLASS_NAME::~DERIVED_CLASS_NAME() {}

const AlgorithmType DERIVED_CLASS_NAME::getType()
{
    return AlgorithmType::BackwardEulerDebugMode;
}

const char* DERIVED_CLASS_NAME::getTypeString()
{
    return "BackwardEulerDebugMode";
}

void DERIVED_CLASS_NAME::Init()
{
    const int allReducesCount = 0;
    hpx_bcast_rsync(neurox::Neuron::SlidingTimeWindow::SetReductionsPerCommStep,
                    &allReducesCount, sizeof(int));
}

void DERIVED_CLASS_NAME::Finalize() {}

double DERIVED_CLASS_NAME::Run()
{
    int commStepSize = Neuron::CommunicationBarrier::commStepSize;
    int totalSteps = Algorithm::getTotalStepsCount();

    hpx_time_t now = hpx_time_now();
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
    double elapsedTime = hpx_time_elapsed_ms(now)/1e3;
    input::Debugger::CompareAllBranches();
    return elapsedTime;
}

void DERIVED_CLASS_NAME::StepBegin(Branch*) {}

void DERIVED_CLASS_NAME::StepEnd(Branch*) {}

void DERIVED_CLASS_NAME::CommStepBegin(Branch*) {}

void DERIVED_CLASS_NAME::CommStepEnd(Branch*) {}

void DERIVED_CLASS_NAME::AfterSpike(Branch*) {}
