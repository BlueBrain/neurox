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

hpx_action_t InitMechanismsAndQuit = 0;
static int InitMechanismsAndQuit_handler()
{
    DebugMessage("neurox::InitMechanismsAndQuit...\n");
    hpx_bcast_rsync(neurox::input::DataLoader::InitMechanisms);
    input::DataLoader::CleanCoreneuronData(false); //clean_ion_global_map=false
    hpx_exit(0,NULL);
}

hpx_action_t Main = 0;
static int Main_handler()
{
    printf("\nneurox::Main (localities: %d, threads/locality: %d, %s)\n",
           hpx_get_num_ranks(), hpx_get_num_threads(), LAYOUT==0 ? "SoA" : "AoS");

    if (!mechanisms) //if not set by neurox::LoadMechsAndQuit
    {
      DebugMessage("neurox::Input::DataLoader::InitMechanisms...\n");
      hpx_bcast_rsync(neurox::input::DataLoader::InitMechanisms);
    }
    DebugMessage("neurox::Input::DataLoader::Init...\n");
    hpx_bcast_rsync(neurox::input::DataLoader::Init);
    DebugMessage("neurox::Input::DataLoader::InitNeurons...\n");
    hpx_bcast_rsync(neurox::input::DataLoader::InitNeurons);
    DebugMessage("neurox::Input::DataLoader::InitNetcons...\n");
    neurox_hpx_call_neurons( neurox::input::DataLoader::InitNetcons);
    DebugMessage("neurox::Input::DataLoader::Finalize...\n");
    hpx_bcast_rsync(neurox::input::DataLoader::Finalize);
    DebugMessage("neurox::Branch::BranchTree::InitLCOs...\n");
    neurox_hpx_call_neurons(Branch::BranchTree::InitLCOs);

    if (neurox::inputParams->outputStatistics)
    {
      DebugMessage("neurox::tools::Statistics::OutputMechanismsDistribution...\n");
      tools::Statistics::OutputMechanismsDistribution();
      DebugMessage("neurox::tools::Statistics::OutputSimulationSize...\n");
      tools::Statistics::OutputSimulationSize();
      //hpx_exit(0,NULL);
    }

    neurox::input::Debugger::CompareMechanismsFunctions();
    neurox::input::Debugger::CompareAllBranches();

    DebugMessage("neurox::Branch::Finitialize...\n");
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
        //TODO for this to work, we have to re-set algorothm in all cpus?
        for (int type = 0; type<AlgorithmType::BenchmarkEnd; type++)
        {
            algorithm = Algorithm::New((AlgorithmType) type);
            algorithm->Init();
            algorithm->PrintStartInfo();
            double timeElapsed = algorithm->Launch();
            totalTimeElapsed += timeElapsed;
            algorithm->Clear();
            delete algorithm;

#ifdef NDEBUG
            //output benchmark info
            printf("csv,%d,%d,%d,%.1f,%.1f,%d,%d,%d,%d,%.2f\n", neurons->size(),
                   hpx_get_num_ranks(), hpx_get_num_threads(), neurons->size() / (double) hpx_get_num_ranks(),
                   inputParams->tstop, algorithm->getType(), inputParams->multiMex ? 1:0,
                   inputParams->branchingDepth, inputParams->allReduceAtLocality ? 1:0, timeElapsed);
            fflush(stdout);
#endif
        }
    }
    else
    {
        algorithm = Algorithm::New(inputParams->algorithm);
        algorithm->Init();
        algorithm->PrintStartInfo();
        totalTimeElapsed = algorithm->Launch();
        algorithm->Clear();
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
        AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::localityNeurons->clear();
        delete AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::localityNeurons;
        AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::localityNeurons = nullptr;
    }

#ifndef NDEBUG
    neurox::input::DataLoader::CleanCoreneuronData();
#endif
    neurox_hpx_unpin;
}

void SetMechanisms(int mechsCount, Mechanism* mechanisms_serial, int * dependenciesIds_serial,
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

void DebugMessage(const char * str)
{
#ifndef NDEBUG
    printf ("%s",str);
    fflush(stdout);
#endif
}

bool ParallelExecution()
{
    return hpx_get_num_ranks()>1;
}

void RegisterHpxActions()
{
    neurox_hpx_register_action(neurox_zero_var_action, neurox::Main);
    neurox_hpx_register_action(neurox_zero_var_action, neurox::InitMechanismsAndQuit);
    neurox_hpx_register_action(neurox_zero_var_action, neurox::Clear);
}

}; //neurox
