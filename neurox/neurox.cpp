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

hpx_action_t Main = 0;
static int Main_handler()
{
    printf("\nneurox::Main (localities: %d, threads/locality: %d, %s)\n",
           hpx_get_num_ranks(), hpx_get_num_threads(), LAYOUT==0 ? "SoA" : "AoS");
    DebugMessage("neurox::Input::DataLoader::Init...\n");
    hpx_bcast_rsync(neurox::input::DataLoader::Init);
    DebugMessage("neurox::Input::DataLoader::InitMechanisms...\n");
    hpx_bcast_rsync(neurox::input::DataLoader::InitMechanisms);
    DebugMessage("neurox::Input::DataLoader::InitNeurons...\n");
    hpx_bcast_rsync(neurox::input::DataLoader::InitNeurons);
    DebugMessage("neurox::Input::DataLoader::InitNetcons...\n");
    neurox_hpx_call_neurons( neurox::input::DataLoader::InitNetcons);
    DebugMessage("neurox::Branch::BranchTree::InitLCOs...\n");
    neurox_hpx_call_neurons(Branch::BranchTree::InitLCOs);
    DebugMessage("neurox::Input::DataLoader::Finalize...\n");
    hpx_bcast_rsync(neurox::input::DataLoader::Finalize);

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

void SetMechanismsDependencies(const char *mechUsed,
                               const int *dependenciesCount, const int * dependenciesIds,
                               const int *successorsCount  , const int * successorsIds)
{
    //make sure mechanisms have already been set
    assert(neurox::mechanismsCount>0 && neurox::mechanisms!= nullptr && neurox::mechanismsMap!=nullptr);

    int successorsIdsOffset=0, dependenciessIdsOffset=0;

    for (int m=0; m<neurox::mechanismsCount; m++)
    {
        Mechanism * mech = mechanisms[m];

        //if this mech is used by this of other locality, mark it as useful
        if (mechUsed[m]) neurox::mechanisms[m]->isUsed = true;

        //if I know about more dependencies dependencies
        if (dependenciesCount!=nullptr)
        {
            assert(dependenciesIds!=nullptr);
            //Less of equal because if its a single dependency (eg on non graph mechs)
            //we are updating it to a new value
            if (mech->dependenciesCount<=dependenciesCount[m])
            {
                delete [] mech->dependencies;
                mech->dependenciesCount = dependenciesCount[m];
                mech->dependencies = new int [dependenciesCount[m]];
                memcpy(mech->dependencies, &dependenciesIds[dependenciessIdsOffset], sizeof(int)*dependenciesCount[m]);
            }
            dependenciessIdsOffset+=dependenciesCount[m];
        }

        //if I know about more successors
        if (successorsCount!=nullptr)
        {
            assert(successorsIds!=nullptr);
            if (mech->successorsCount<=successorsCount[m])
            {
                delete [] mech->successors;
                mech->successorsCount = successorsCount[m];
                mech->successors = new int [successorsCount[m]];
                memcpy(mech->successors, &successorsIds[successorsIdsOffset], sizeof(int)*successorsCount[m]);
            }
            successorsIdsOffset+=successorsCount[m];
        }
    }

    //initializes parent ion index
    for (int m=0; m<mechanismsCount; m++)
    {
      Mechanism * mech = mechanisms[m];
      if (inputParams->multiMex)
      {
        for (int d=0; d<mech->dependenciesCount; d++)
        {
          Mechanism * parent = GetMechanismFromType(mech->dependencies[d]);
          if (strcmp("SK_E2", mech->membFunc.sym)==0 && strcmp("ca_ion", parent->membFunc.sym)==0) continue; //TODO hard coded exception
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
    neurox_hpx_register_action(neurox_zero_var_action, neurox::Clear);
}

}; //neurox
