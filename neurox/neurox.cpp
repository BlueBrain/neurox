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
neurox::algorithms::Algorithm * algorithm = nullptr;

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

    double totalTimeElapsed = 0;
    if (inputParams->algorithm == AlgorithmType::All)
    {
        for (int type = 0; type<AlgorithmType::BenchmarkEnd; type++)
        {
            algorithm = Algorithm::New((AlgorithmType) type);
            algorithm->Init();
            algorithm->PrintStartInfo();
            double timeElapsed = algorithm->Run();
            totalTimeElapsed += timeElapsed;
            algorithm->Finalize();
            delete algorithm;

#ifdef NDEBUG
    //output benchmark info
    double timeElapsed = hpx_time_elapsed_ms(now)/1e3;
    printf("csv,%d,%d,%d,%.1f,%.1f,%d,%d,%d,%d,%.2f\n", neurons->size(), hpx_get_num_ranks(),
        hpx_get_num_threads(), neurons->size() / (double) hpx_get_num_ranks(), inputParams->tstop,
        algorithm->getType(), inputParams->multiMex ? 1:0, inputParams->branchingDepth,
        inputParams->allReduceAtLocality ? 1:0, timeElapsed);
    fflush(stdout);
#endif
        }
    }
    else
    {
        algorithm = Algorithm::New(inputParams->algorithm);
        algorithm->Init();
        algorithm->PrintStartInfo();
        totalTimeElapsed = algorithm->Run();
        algorithm->Finalize();
        delete algorithm;
    }

    printf("neurox::end (%d neurons, biological time: %.3f secs, solver time: %.3f secs).\n",
           neurons->size(), inputParams->tstop/1000.0, totalTimeElapsed);

    neurox_hpx_call_neurons(Branch::Clear);
    hpx_bcast_rsync(neurox::Clear);
    hpx_exit(0,NULL);
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
    neurox_hpx_register_action(neurox_several_vars_action, neurox::SetMechanisms);
    neurox_hpx_register_action(neurox_several_vars_action, neurox::SetMechanismsGlobalVars);
}

}; //neurox
