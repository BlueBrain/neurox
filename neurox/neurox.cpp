#include "neurox/neurox.h"
#include <cstring>
#include <map>
#include <fstream>
#include <iostream>

#include "coreneuron/nrniv/nrniv_decl.h"
#include "coreneuron/nrniv/nrn_stats.h"

using namespace neurox::algorithms;

namespace neurox
{

std::vector<hpx_t> * neurons = nullptr;
int mechanismsCount = -1;
int * mechanismsMap = nullptr;
neurox::Mechanism ** mechanisms = nullptr;
neurox::tools::CmdLineParser * inputParams = nullptr;

Mechanism * GetMechanismFromType(int type) {
    assert(mechanismsMap[type]!=-1);
    return mechanisms[mechanismsMap[type]];
}

void SetMechanisms2(int mechsCount, Mechanism* mechanisms_serial, int * dependenciesIds_serial,
                    int * successorsIds_serial, char * sym_serial)
{
    neurox::mechanismsCount = mechsCount;
    neurox::mechanisms = new Mechanism*[mechsCount];
    int offsetSuccessors=0, offsetDependencies=0;
    int offsetSym=0;
    int maxMechType=-1;
    for (int m=0; m<mechanismsCount; m++)
    {
        Mechanism & mech = mechanisms_serial[m];
        int * dependenciesIds = mech.dependenciesCount == 0 ? nullptr : &dependenciesIds_serial[offsetDependencies];
        int * successorsIds = mech.successorsCount == 0 ? nullptr : &successorsIds_serial[offsetSuccessors];
        char * sym = mech.symLength == 0 ? nullptr : &sym_serial[offsetSym];
        mechanisms[m] = new Mechanism(
                    mech.type, mech.dataSize, mech.pdataSize,
                    mech.isArtificial, mech.pntMap, mech.isIon,
                    mech.symLength, sym,
                    mech.dependenciesCount, dependenciesIds,
                    mech.successorsCount, successorsIds);
        offsetSuccessors +=  mech.successorsCount;
        offsetDependencies +=  mech.dependenciesCount;
        offsetSym += mech.symLength;
        if (mech.type > maxMechType)
            maxMechType = mech.type;
    }

    //initializes map of mechanisms ids to offset
    neurox::mechanismsMap = new int[maxMechType+1];
    for (int i=0; i<maxMechType+1; i++)
        neurox::mechanismsMap[i]=-1;
    for (int m=0; m<mechanismsCount; m++)
        neurox::mechanismsMap[mechanisms[m]->type]=m;

    //initializes parent ion index
    for (int m=0; m<mechanismsCount; m++)
    {
      Mechanism * mech = mechanisms[m];
      mech->dependencyIonIndex = Mechanism::Ion::no_ion;
      if (inputParams->multiMex)
      {
        for (int d=0; d<mech->dependenciesCount; d++)
        {
          Mechanism * parent = GetMechanismFromType(mech->dependencies[d]);
          if (strcmp("SK_E2", mech->sym)==0 && strcmp("ca_ion", parent->sym)==0) continue; //TODO hard coded exception
          if (parent->GetIonIndex() < Mechanism::Ion::size_writeable_ions)
              mech->dependencyIonIndex = parent->GetIonIndex();
        }
      }
    }
}

hpx_action_t SetAlgorithmVariables = 0;
int SetAlgorithmVariables_handler(const AlgorithmType * algorithm_ptr, const size_t)
{
    neurox_hpx_pin(uint64_t);
    inputParams->algorithm = *algorithm_ptr;
    Neuron::SlidingTimeWindow::reductionsPerCommStep = 0;
    if (*algorithm_ptr==AlgorithmType::BackwardEulerSlidingTimeWindow)
        Neuron::SlidingTimeWindow::reductionsPerCommStep = 2;
    else if ( *algorithm_ptr==AlgorithmType::BackwardEulerAllReduce)
        Neuron::SlidingTimeWindow::reductionsPerCommStep = 1;
    neurox_hpx_unpin;
}

hpx_action_t SetMechanisms = 0;
int SetMechanisms_handler(const int nargs, const void *args[], const size_t sizes[])
{
    /**
     * nargs=4 where:
     * args[0] = array of all mechanisms info
     * args[1] = array of all mechanisms dependencies (parents in mechanisms tree)
     * args[2] = array of all mechanisms successors (children in mechanisms tree)
     * args[3] = array of all mechanisms names (sym)
     */
    neurox_hpx_pin(uint64_t);
    assert(nargs==4);
    int mechanismsCount = sizes[0]/sizeof(Mechanism);
    SetMechanisms2(mechanismsCount, (Mechanism*) args[0], (int*) args[1], (int*) args[2], (char*) args[3]);
    neurox_hpx_unpin;
}

hpx_action_t SetMechanismsGlobalVars = 0;
int SetMechanismsGlobalVars_handler(const int nargs, const void *args[], const size_t sizes[])
{
    /**
     * nargs=3 where:
     * args[0] = celsius
     * args[1] = nrn_ion_global_map_size
     * args[2] = flag if types has entry in nrn_ion_global_map, 1 or 0 (ie not null)
     * args[3] = nrn_ion_global_map
     */

    neurox_hpx_pin(uint64_t);
    assert(nargs==4);
    double new_celsius = *(double*) args[0];
    int mechsCount = *(int*) args[1];
    unsigned char * mechHasEntryInIonMap = (unsigned char*) args[2];
    double * ionGlobalMapInfo = (double*) args[3];

    int mechsOffset=0;
    if (nrn_ion_global_map!=NULL) //this machine has the data, compare
    {
        assert (celsius == new_celsius);
        assert (mechsCount == nrn_ion_global_map_size);
        for (int i=0; i<mechsCount; i++)
        {
            if (mechHasEntryInIonMap[i]==0)
            {
                assert (nrn_ion_global_map[i] == NULL);
            }
            else
            {
                assert (mechHasEntryInIonMap[i]==1);
                assert (nrn_ion_global_map[i][0] = ionGlobalMapInfo[mechsOffset+0]);
                assert (nrn_ion_global_map[i][1] = ionGlobalMapInfo[mechsOffset+1]);
                assert (nrn_ion_global_map[i][2] = ionGlobalMapInfo[mechsOffset+2]);
                mechsOffset += 3;
            }
        }
    }
    else //this machine does not have the data, create it
    {
        celsius = new_celsius;
        nrn_ion_global_map = new double*[mechsCount];
        for (int i=0; i<mechsCount; i++)
        {
            if (!mechHasEntryInIonMap[i])
            {
                nrn_ion_global_map[i] = NULL;
            }
            else
            {
                nrn_ion_global_map[i] = new double[3];
                nrn_ion_global_map[i][0] = ionGlobalMapInfo[mechsOffset+0];
                nrn_ion_global_map[i][1] = ionGlobalMapInfo[mechsOffset+1];
                nrn_ion_global_map[i][2] = ionGlobalMapInfo[mechsOffset+2];
                mechsOffset+=3;
            }
        }
    }
    neurox_hpx_unpin;
}

void DebugMessage(const char * str)
{
#ifndef NDEBUG
    printf ("%s",str);
    fflush(stdout);
#endif
}


hpx_action_t Main = 0;
static int Main_handler()
{
    printf("\nneurox::main (localities: %d, threads/locality: %d, %s)\n",
           hpx_get_num_ranks(), hpx_get_num_threads(), LAYOUT==0 ? "SoA" : "AoS");
    if (hpx_get_num_ranks()>1 && !inputParams->parallelDataLoading)
    {
      DebugMessage("ERROR: add the -m or --mpi argument for parallel data loading\n");
      hpx_exit(0, NULL);
    }
    DebugMessage("neurox::Input::DataLoader::init...\n");
    hpx_bcast_rsync(neurox::input::DataLoader::Init);
    DebugMessage("neurox::Input::DataLoader::initMechanisms...\n");
    hpx_bcast_rsync(neurox::input::DataLoader::InitMechanisms);
    DebugMessage("neurox::Input::DataLoader::initNeurons...\n");
    hpx_bcast_rsync(neurox::input::DataLoader::InitNeurons);
    DebugMessage("neurox::Input::DataLoader::initNetcons...\n");
    neurox_hpx_call_neurons( neurox::input::DataLoader::InitNetcons, nullptr, 0);
    DebugMessage("neurox::Input::DataLoader::finalize...\n");
    hpx_bcast_rsync(neurox::input::DataLoader::Finalize);
    DebugMessage("neurox::Branch::BranchTree::initLCOs...\n");
    neurox_hpx_call_neurons(Branch::BranchTree::InitLCOs);

    if (neurox::inputParams->outputStatistics)
    {
      DebugMessage("neurox::Misc::Statistics::outputMechanismsDistribution...\n");
      tools::Statistics::OutputMechanismsDistribution();
      DebugMessage("neurox::Misc::Statistics::outputSimulationSize...\n");
      tools::Statistics::OutputSimulationSize();
      //hpx_exit(0,NULL);
    }

    neurox::input::Debugger::CompareMechanismsFunctionPointers();
    neurox::input::Debugger::CompareAllBranches();

    DebugMessage("neurox::Branch::finitialize...\n");
    neurox_hpx_call_neurons(Branch::Finitialize);
#ifndef NDEBUG
    hpx_bcast_rsync(neurox::input::Debugger::Finitialize);
    neurox::input::Debugger::CompareAllBranches();
#endif

    DebugMessage("neurox::Branch::threadTableCheck...\n");
    neurox_hpx_call_neurons(Branch::ThreadTableCheck);
#ifndef NDEBUG
    hpx_bcast_rsync(neurox::input::Debugger::ThreadTableCheck);
    neurox::input::Debugger::CompareAllBranches();
#endif

    //subscribe to the all-reduce LCOs
    static hpx_t * allreduces = nullptr;
    static int maxReductionsPerCommStep = 0;
    if (inputParams->algorithm == AlgorithmType::All
     || inputParams->algorithm == AlgorithmType::BackwardEulerSlidingTimeWindow
     || inputParams->algorithm == AlgorithmType::BackwardEulerAllReduce)
    {
        DebugMessage("neurox::Neuron::SlidingTimeWindow::init...\n");
        maxReductionsPerCommStep = inputParams->algorithm == AlgorithmType::BackwardEulerAllReduce ? 1 : 2;
        hpx_bcast_rsync(neurox::Neuron::SlidingTimeWindow::SetReductionsPerCommStep, &maxReductionsPerCommStep, sizeof(int));

        allreduces = new hpx_t[maxReductionsPerCommStep];
        for (int i=0; i<maxReductionsPerCommStep; i++)
            allreduces[i] = hpx_process_collective_allreduce_new(0, Neuron::SlidingTimeWindow::Init, Neuron::SlidingTimeWindow::Reduce);

        if (inputParams->allReduceAtLocality)
            hpx_bcast_rsync(neurox::Neuron::SlidingTimeWindow::AllReduceLocality::SubscribeAllReduce,
                            allreduces, sizeof(hpx_t)*maxReductionsPerCommStep);
        else
            neurox_hpx_call_neurons(Neuron::SlidingTimeWindow::SubscribeAllReduce,
                            allreduces, sizeof(hpx_t)*maxReductionsPerCommStep);

        for (int i=0; i<maxReductionsPerCommStep; i++)
            hpx_process_collective_allreduce_subscribe_finalize(allreduces[i]);
    }

    hpx_time_t now = hpx_time_now();
    if (inputParams->algorithm == AlgorithmType::All)
    {
        RunAlgorithm(AlgorithmType::BackwardEulerAllReduce );
        RunAlgorithm(AlgorithmType::BackwardEulerSlidingTimeWindow);
        RunAlgorithm(AlgorithmType::BackwardEulerTimeDependencyLCO);
    }
    else
        RunAlgorithm(inputParams->algorithm);

    printf("neurox::end (%d neurons, biological time: %.3f secs, solver time: %.3f secs).\n",
           neurons->size(), inputParams->tstop/1000.0, hpx_time_elapsed_ms(now)/1e3);

    //Clean up all-reduce LCOs
    if ( inputParams->algorithm == All
      || inputParams->algorithm == AlgorithmType::BackwardEulerSlidingTimeWindow
      || inputParams->algorithm == AlgorithmType::BackwardEulerAllReduce)
    {
        if (inputParams->allReduceAtLocality)
            hpx_bcast_rsync(neurox::Neuron::SlidingTimeWindow::AllReduceLocality::UnsubscribeAllReduce,
                            allreduces, sizeof(hpx_t)*maxReductionsPerCommStep);
        else
            neurox_hpx_call_neurons(Neuron::SlidingTimeWindow::UnsubscribeAllReduce, allreduces, sizeof(hpx_t)*maxReductionsPerCommStep);

        for (int i=0; i<maxReductionsPerCommStep; i++)
            hpx_process_collective_allreduce_delete(allreduces[i]);
        delete [] allreduces; allreduces=nullptr;
    }

    neurox_hpx_call_neurons(Branch::Clear);
    hpx_bcast_rsync(neurox::Clear);
    hpx_exit(0,NULL);
}

void RunAlgorithm(AlgorithmType algorithm)
{
    int totalSteps = (inputParams->tstop - inputParams->tstart) / inputParams->dt;
    printf("neurox::Algorithm::%s (%d neurons, t=%.03f secs, dt=%.03f milisecs, %d steps)\n",
            algorithm == AlgorithmType::BackwardEulerDebugMode  ? "BackwardEulerDebugWithCommBarrier" :
           (algorithm == AlgorithmType::All                                ? "ALL" :
           (algorithm == AlgorithmType::BackwardEulerAllReduce  ? "BackwardEulerWithAllReduceBarrier" :
           (algorithm == AlgorithmType::BackwardEulerSlidingTimeWindow ? "BackwardEulerWithSlidingTimeWindow" :
           (algorithm == AlgorithmType::BackwardEulerTimeDependencyLCO ? "BackwardEulerWithTimeDependencyLCO" :
            "unknown" )))), neurons->size(), inputParams->tstop/1000, inputParams->dt, totalSteps);
    fflush(stdout);

    AlgorithmType previousAlgorithm = inputParams->algorithm;
    hpx_bcast_rsync(neurox::SetAlgorithmVariables, &algorithm, sizeof(AlgorithmType));

    hpx_t mainLCO = hpx_lco_and_new(neurons->size());

#ifdef NDEBUG //benchmark info
    hpx_time_t now = hpx_time_now();
#endif

    if (inputParams->algorithm == AlgorithmType::BackwardEulerDebugMode)
    {
        int commStepSize = Neuron::CommunicationBarrier::commStepSize;
        for (int s=0; s<totalSteps; s+=Neuron::CommunicationBarrier::commStepSize)
        {
            #ifdef NEUROX_TIME_STEPPING_VERBOSE
              if (hpx_get_my_rank()==0)
                DebugMessage(std::string("-- t="+std::to_string(inputParams->dt*s)+" ms\n").c_str());
            #endif

            //Reduction at locality is not implemented (this mode is for debugging only)
            neurox_hpx_call_neurons_lco(Branch::BackwardEuler, mainLCO, &commStepSize, sizeof(int));

            #ifndef NDEBUG
              if (inputParams->parallelDataLoading) //if parallel execution... spike exchange
                hpx_bcast_rsync(neurox::input::Debugger::NrnSpikeExchange);
            #endif
        }
    }
    else
    {
        if (inputParams->allReduceAtLocality &&
            (  inputParams->algorithm == AlgorithmType::BackwardEulerAllReduce
            || inputParams->algorithm == AlgorithmType::BackwardEulerSlidingTimeWindow))
            hpx_bcast_rsync(Branch::BackwardEulerOnLocality, &totalSteps, sizeof(int));
        else
            neurox_hpx_call_neurons_lco(Branch::BackwardEuler, mainLCO, &totalSteps, sizeof(int));
    }

#ifdef NDEBUG
    //output benchmark info
    double timeElapsed = hpx_time_elapsed_ms(now)/1e3;
    printf("csv,%d,%d,%d,%.1f,%.1f,%d,%d,%d,%d,%.2f\n", neurons->size(), hpx_get_num_ranks(),
        hpx_get_num_threads(), neurons->size() / (double) hpx_get_num_ranks(), inputParams->tstop,
        inputParams->algorithm, inputParams->multiMex ? 1:0, inputParams->branchingDepth,
        inputParams->allReduceAtLocality ? 1:0, timeElapsed);
    fflush(stdout);
#else
    //compare final results
    if (inputParams->branchingDepth==0)
    if (inputParams->algorithm != AlgorithmType::BackwardEulerDebugMode //not fixed comm barrier
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
    hpx_lco_delete_sync(mainLCO);
}

hpx_action_t Clear = 0;
int Clear_handler()
{
    neurox_hpx_pin(uint64_t);
    delete [] mechanisms;
    if (neurox::neurons)
    {
        neurox::neurons->clear();
        delete neurox::neurons;
        neurox::neurons = nullptr;
    }

    if (inputParams->allReduceAtLocality)
    {
        Neuron::SlidingTimeWindow::AllReduceLocality::localityNeurons->clear();
        delete Neuron::SlidingTimeWindow::AllReduceLocality::localityNeurons;
        Neuron::SlidingTimeWindow::AllReduceLocality::localityNeurons = nullptr;
    }

#ifndef NDEBUG
    neurox::input::DataLoader::CleanCoreneuronData();
#endif
    neurox_hpx_unpin;
}

void RegisterHpxActions()
{
    neurox_hpx_register_action(neurox_zero_var_action,     neurox::Main);
    neurox_hpx_register_action(neurox_zero_var_action,     neurox::Clear);
    neurox_hpx_register_action(neurox_single_var_action,   neurox::SetAlgorithmVariables);
    neurox_hpx_register_action(neurox_several_vars_action, neurox::SetMechanisms);
    neurox_hpx_register_action(neurox_several_vars_action, neurox::SetMechanismsGlobalVars);
}

}; //neurox
