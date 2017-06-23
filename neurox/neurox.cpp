#include "neurox/neurox.h"
#include <cstring>
#include <map>
#include <fstream>
#include <iostream>

#include "coreneuron/nrniv/nrniv_decl.h"
#include "coreneuron/nrniv/nrn_stats.h"

namespace neurox
{

std::vector<hpx_t> * neurons = nullptr;
int mechanismsCount = -1;
int * mechanismsMap = nullptr;
neurox::Mechanism ** mechanisms = nullptr;
Input::InputParams * inputParams = nullptr;

Mechanism * getMechanismFromType(int type) {
    assert(mechanismsMap[type]!=-1);
    return mechanisms[mechanismsMap[type]];
}

void setMechanisms2(int mechsCount, Mechanism* mechanisms_serial, int * dependenciesIds_serial,
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
          Mechanism * parent = getMechanismFromType(mech->dependencies[d]);
          if (strcmp("SK_E2", mech->sym)==0 && strcmp("ca_ion", parent->sym)==0) continue; //TODO hard coded exception
          if (parent->getIonIndex() < Mechanism::Ion::size_writeable_ions)
              mech->dependencyIonIndex = parent->getIonIndex();
        }
      }
    }
}

hpx_action_t setAlgorithmVariables = 0;
int setAlgorithmVariables_handler(const Algorithm * algorithm_ptr, const size_t)
{
    neurox_hpx_pin(uint64_t);
    inputParams->algorithm = *algorithm_ptr;
    Neuron::SlidingTimeWindow::reductionsPerCommStep = 0;
    if (*algorithm_ptr==Algorithm::BackwardEulerWithSlidingTimeWindow)
        Neuron::SlidingTimeWindow::reductionsPerCommStep = 2;
    else if ( *algorithm_ptr==Algorithm::BackwardEulerWithAllReduceBarrier)
        Neuron::SlidingTimeWindow::reductionsPerCommStep = 1;
    neurox_hpx_unpin;
}

hpx_action_t setMechanisms = 0;
int setMechanisms_handler(const int nargs, const void *args[], const size_t sizes[])
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
    setMechanisms2(mechanismsCount, (Mechanism*) args[0], (int*) args[1], (int*) args[2], (char*) args[3]);
    neurox_hpx_unpin;
}

hpx_action_t setMechanismsGlobalVars = 0;
int setMechanismsGlobalVars_handler(const int nargs, const void *args[], const size_t sizes[])
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

void message(const char * str)
{
#ifndef NDEBUG
    printf ("%s",str);
    fflush(stdout);
#endif
}


hpx_action_t main = 0;
static int main_handler()
{
    printf("\nneurox::main (localities: %d, threads/locality: %d)\n", hpx_get_num_ranks(), hpx_get_num_threads());
    if (hpx_get_num_ranks()>1 && !inputParams->parallelDataLoading)
    {
      message("ERROR: add the -m or --mpi argument for parallel data loading\n");
      hpx_exit(0, NULL);
    }
    message("neurox::Input::Coreneuron::DataLoader::init...\n");
    hpx_bcast_rsync(neurox::Input::Coreneuron::DataLoader::init);
    message("neurox::Input::Coreneuron::DataLoader::initMechanisms...\n");
    hpx_bcast_rsync(neurox::Input::Coreneuron::DataLoader::initMechanisms);
    message("neurox::Input::Coreneuron::DataLoader::initNeurons...\n");
    hpx_bcast_rsync(neurox::Input::Coreneuron::DataLoader::initNeurons);
    message("neurox::Input::Coreneuron::DataLoader::initNetcons...\n");
    neurox_hpx_call_neurons( neurox::Input::Coreneuron::DataLoader::initNetcons);
    message("neurox::Input::Coreneuron::DataLoader::finalize...\n");
    hpx_bcast_rsync(neurox::Input::Coreneuron::DataLoader::finalize);
    message("neurox::Branch::BranchTree::initLCOs...\n");
    neurox_hpx_call_neurons(Branch::BranchTree::initLCOs);

    if (neurox::inputParams->outputStatistics)
    {
      message("neurox::Misc::Statistics::printMechanismsDistribution...\n");
      Misc::Statistics::printMechanismsDistribution();
      message("neurox::Misc::Statistics::printSimulationSize...\n");
      Misc::Statistics::printSimulationSize();
      //hpx_exit(0,NULL);
    }

    neurox::Input::Coreneuron::Debugger::compareMechanismsFunctionPointers();
    neurox::Input::Coreneuron::Debugger::compareAllBranches();

    message("neurox::Branch::finitialize...\n");
    neurox_hpx_call_neurons(Branch::finitialize);
#ifndef NDEBUG
    hpx_bcast_rsync(neurox::Input::Coreneuron::Debugger::finitialize);
    neurox::Input::Coreneuron::Debugger::compareAllBranches();
#endif

    message("neurox::Branch::threadTableCheck...\n");
    neurox_hpx_call_neurons(Branch::threadTableCheck);
#ifndef NDEBUG
    hpx_bcast_rsync(neurox::Input::Coreneuron::Debugger::threadTableCheck);
    neurox::Input::Coreneuron::Debugger::compareAllBranches();
#endif

    //subscribe to the all-reduce LCOs
    static hpx_t * allreduces = nullptr;
    static int maxReductionsPerCommStep = 0;
    if (inputParams->algorithm == Algorithm::ALL
     || inputParams->algorithm == Algorithm::BackwardEulerWithSlidingTimeWindow
     || inputParams->algorithm == Algorithm::BackwardEulerWithAllReduceBarrier)
    {
        message("neurox::Neuron::SlidingTimeWindow::init...\n");
        maxReductionsPerCommStep = inputParams->algorithm == Algorithm::BackwardEulerWithAllReduceBarrier ? 1 : 2;
        hpx_bcast_rsync(neurox::Neuron::SlidingTimeWindow::setReductionsPerCommStep, &maxReductionsPerCommStep, sizeof(int));

        allreduces = new hpx_t[maxReductionsPerCommStep];
        for (int i=0; i<maxReductionsPerCommStep; i++)
            allreduces[i] = hpx_process_collective_allreduce_new(0, Neuron::SlidingTimeWindow::init, Neuron::SlidingTimeWindow::reduce);

        if (inputParams->allReduceAtLocality)
            hpx_bcast_rsync(neurox::Neuron::SlidingTimeWindow::AllReduceLocality::subscribeAllReduce,
                            allreduces, sizeof(hpx_t)*maxReductionsPerCommStep);
        else
            neurox_hpx_call_neurons(Neuron::SlidingTimeWindow::subscribeAllReduce,
                            allreduces, sizeof(hpx_t)*maxReductionsPerCommStep);

        for (int i=0; i<maxReductionsPerCommStep; i++)
            hpx_process_collective_allreduce_subscribe_finalize(allreduces[i]);
    }

    int totalSteps = (inputParams->tstop - inputParams->tstart) / inputParams->dt;
    printf("neurox::Algorithm::%s (%d neurons, t=%.03f secs, dt=%.03f milisecs, %d steps)\n",
            inputParams->algorithm == Algorithm::BackwardEulerDebugWithCommBarrier  ? "BackwardEulerDebugWithCommBarrier" :
           (inputParams->algorithm == Algorithm::ALL                                ? "ALL" :
           (inputParams->algorithm == Algorithm::BackwardEulerWithAllReduceBarrier  ? "BackwardEulerWithAllReduceBarrier" :
           (inputParams->algorithm == Algorithm::BackwardEulerWithSlidingTimeWindow ? "BackwardEulerWithSlidingTimeWindow" :
           (inputParams->algorithm == Algorithm::BackwardEulerWithTimeDependencyLCO ? "BackwardEulerWithTimeDependencyLCO" :
            "UNKNOWN" )))), neurons->size(), inputParams->tstop/1000, inputParams->dt, totalSteps);
    fflush(stdout);

    hpx_time_t now = hpx_time_now();
    if (inputParams->algorithm == Algorithm::ALL)
    {
        //TODO Why running all three is slower (for BackwardEulerWithAllReduceBarrier) than running individually?
        runAlgorithm(Algorithm::BackwardEulerWithAllReduceBarrier );
        runAlgorithm(Algorithm::BackwardEulerWithSlidingTimeWindow);
        runAlgorithm(Algorithm::BackwardEulerWithTimeDependencyLCO);
    }
    else
        runAlgorithm(inputParams->algorithm);

    printf("neurox::end (%d neurons, biological time: %.3f secs, solver time: %.3f secs).\n",
           neurons->size(), inputParams->tstop/1000.0, hpx_time_elapsed_ms(now)/1e3);

    //Clean up all-reduce LCOs
    if ( inputParams->algorithm == ALL
      || inputParams->algorithm == Algorithm::BackwardEulerWithSlidingTimeWindow
      || inputParams->algorithm == Algorithm::BackwardEulerWithAllReduceBarrier)
    {
        if (inputParams->allReduceAtLocality)
            hpx_bcast_rsync(neurox::Neuron::SlidingTimeWindow::AllReduceLocality::unsubscribeAllReduce,
                            allreduces, sizeof(hpx_t)*maxReductionsPerCommStep);
        else
            neurox_hpx_call_neurons(Neuron::SlidingTimeWindow::unsubscribeAllReduce, allreduces, sizeof(hpx_t)*maxReductionsPerCommStep);

        for (int i=0; i<maxReductionsPerCommStep; i++)
            hpx_process_collective_allreduce_delete(allreduces[i]);
        delete [] allreduces; allreduces=nullptr;
    }

    neurox_hpx_call_neurons(Branch::clear);
    hpx_bcast_rsync(neurox::clear);
    hpx_exit(0,NULL);
}

void runAlgorithm(Algorithm algorithm)
{
    Algorithm previousAlgorithm = inputParams->algorithm;
    hpx_bcast_rsync(neurox::setAlgorithmVariables, &algorithm, sizeof(Algorithm));

    int totalSteps = (inputParams->tstop - inputParams->tstart) / inputParams->dt;
    hpx_t mainLCO = hpx_lco_and_new(neurons->size());

#ifdef NDEBUG //benchmark info
    hpx_time_t now = hpx_time_now();
#endif

    size_t neuronsCount = neurons->size();
    if (inputParams->algorithm == Algorithm::BackwardEulerDebugWithCommBarrier)
    {
        int commStepSize = Neuron::CommunicationBarrier::commStepSize;
        for (int s=0; s<totalSteps; s+=Neuron::CommunicationBarrier::commStepSize)
        {
            #ifdef NEUROX_TIME_STEPPING_VERBOSE
              if (hpx_get_my_rank()==0)
                message(std::string("-- t="+std::to_string(inputParams->dt*s)+" ms\n").c_str());
            #endif

            //Reduction at locality is not implemented (this mode is for debugging only)
            neurox_hpx_call_neurons_lco(Branch::backwardEuler, mainLCO, &commStepSize, sizeof(int));

            #ifndef NDEBUG
              if (inputParams->parallelDataLoading) //if parallel execution... spike exchange
                hpx_bcast_rsync(neurox::Input::Coreneuron::Debugger::nrnSpikeExchange);
            #endif
        }
    }
    else
    {
        if (inputParams->allReduceAtLocality &&
            (  inputParams->algorithm == Algorithm::BackwardEulerWithAllReduceBarrier
            || inputParams->algorithm == Algorithm::BackwardEulerWithSlidingTimeWindow))
            hpx_bcast_rsync(Branch::backwardEulerOnLocality, &totalSteps, sizeof(int));
        else
            neurox_hpx_call_neurons_lco(Branch::backwardEuler, mainLCO, &totalSteps, sizeof(int));
    }

#ifdef NDEBUG
    //output benchmark info
    double elapsed = hpx_time_elapsed_ms(now)/1e3;
    printf("csv,%d,%d,%d,%.1f,%.1f,%d,%d,%d,%d,%.2f\n", neurons->size(), hpx_get_num_ranks(),
        hpx_get_num_threads(), neurons->size() / (double) hpx_get_num_ranks(), inputParams->tstop,
        inputParams->algorithm, inputParams->multiMex ? 1:0, inputParams->multiSplix ? 1:0,
        inputParams->allReduceAtLocality ? 1:0, elapsed);
    fflush(stdout);
#else
    //compare final results
    if (!inputParams->algorithm == Algorithm::BackwardEulerDebugWithCommBarrier //not fixed comm barrier
    && inputParams->parallelDataLoading) //and not serial
    {
        //re-run whole simulation and comparae final result
        message("neurox::re-running simulation in Coreneuron to compare final result...\n");
        fflush(stdout);
        int commStepSize = Neuron::CommunicationBarrier::commStepSize;
        for (int s=0; s<totalSteps; s+=Neuron::CommunicationBarrier::commStepSize)
        {
            hpx_bcast_rsync(neurox::Input::Coreneuron::Debugger::fixedStepMinimal, &commStepSize, sizeof(int));
            hpx_bcast_rsync(neurox::Input::Coreneuron::Debugger::nrnSpikeExchange);
        }
    }
    neurox::Input::Coreneuron::Debugger::compareAllBranches();
#endif

    hpx_bcast_rsync(neurox::setAlgorithmVariables, &previousAlgorithm, sizeof(Algorithm));
    hpx_lco_delete_sync(mainLCO);
}

hpx_action_t clear = 0;
int clear_handler()
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
    neurox::Input::Coreneuron::DataLoader::cleanCoreneuronData();
#endif
    neurox_hpx_unpin;
}

void registerHpxActions()
{
    neurox_hpx_register_action(0,neurox::main);
    neurox_hpx_register_action(0,neurox::clear);
    neurox_hpx_register_action(1,neurox::setAlgorithmVariables);
    neurox_hpx_register_action(2,neurox::setMechanisms);
    neurox_hpx_register_action(2,neurox::setMechanismsGlobalVars);
}

}; //namespace neurox
