#include "neurox/neurox.h"
#include <cstring>
#include <map>
#include <fstream>

#include "nrniv/nrniv_decl.h"
#include "nrniv/nrn_stats.h"

#ifndef NEUROX_LOGO_PATH
  #define ""
#endif

namespace neurox
{

std::vector<hpx_t> * neurons   = nullptr;
std::vector<hpx_t> * myNeurons = nullptr;
int mechanismsCount=-1;
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
      mech->dependencyIonIndex = Branch::MechanismsGraph::IonIndex::no_ion;
      if (inputParams->multiMex)
      {
        for (int d=0; d<mech->dependenciesCount; d++)
        {
          Mechanism * parent = getMechanismFromType(mech->dependencies[d]);
          if (strcmp("SK_E2", mech->sym)==0 && strcmp("ca_ion", parent->sym)==0) continue; //TODO hard coded exception
          mech->dependencyIonIndex = parent->getIonIndex();
        }
      }
    }
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

static void displayLogo()
{
    string getcontent;
    ifstream openfile (NEUROX_LOGO_PATH);
    if(openfile.is_open())
        while(getline(openfile, getcontent))
            std::cout << getcontent << endl;
    else printf("Warning: can't open logo file %s. Ignoring..\n", NEUROX_LOGO_PATH);
    openfile.close();
}

void message(char * str)
{
#ifdef  NEUROX_TIME_STEPPING_VERBOSE
    printf ("%s",str);
    fflush(stdout);
#endif
}

void compareAllBranches(int neuronsCount)
{
#if !defined(NDEBUG) && defined(CORENEURON_H)
    message("neurox::Input::CoreNeuron::Debugger::compareBranch...\n");
    hpx_par_for_sync( [&] (int i, void*) -> int
    {  return hpx_call_sync(neurox::neurons->at(i), Input::Coreneuron::Debugger::compareBranch, HPX_NULL, 0);
    }, 0, neuronsCount, NULL);
#endif
}

hpx_action_t main = 0;
static int main_handler()
{
#ifdef NDEBUG
    neurox::displayLogo();
#endif
    printf("neurox::main (localities: %d, threads/locality: %d)\n", hpx_get_num_ranks(), hpx_get_num_threads());
    fflush(stdout);
    message("neurox::Input::Coreneuron::DataLoader::init...\n");
    hpx_bcast_rsync(neurox::Input::Coreneuron::DataLoader::init);
    message("neurox::Input::Coreneuron::DataLoader::initMechanisms...\n");
    hpx_bcast_rsync(neurox::Input::Coreneuron::DataLoader::initMechanisms);
    message("neurox::Input::Coreneuron::DataLoader::initNeurons...\n");
    hpx_bcast_rsync(neurox::Input::Coreneuron::DataLoader::initNeurons);
    message("neurox::Input::Coreneuron::DataLoader::initNetcons...\n");
    hpx_par_for_sync( [&] (int i, void*) -> int
    {  return hpx_call_sync(neurox::neurons->at(i), neurox::Input::Coreneuron::DataLoader::initNetcons, HPX_NULL, 0);
    }, 0, neurox::neurons->size(), NULL);
    message("neurox::Input::Coreneuron::DataLoader::clean...\n");
    hpx_bcast_rsync(neurox::Input::Coreneuron::DataLoader::clean);

    int neuronsCount = neurons->size();
    if (neurox::inputParams->outputStatistics)
    {
      message("neurox::Misc::Statistics::printMechanismsDistribution...\n");
      Misc::Statistics::printMechanismsDistribution();
      message("neurox::Misc::Statistics::printSimulationSize...\n");
      Misc::Statistics::printSimulationSize();
      //hpx_exit(0,NULL);
    }

    if (inputParams->algorithm == Algorithm::BackwardEulerSlidingTimeWindow)
    {
        //hpx_t slidingWindowAllReduce = hpx_process_collective_allreduce_new(sizeof(int), _reset, _op);
    }

    message("neurox::Branch::BranchTree::initLCOs...\n");
    hpx_par_for_sync( [&] (int i, void*) -> int
    {  return hpx_call_sync(neurox::neurons->at(i), Branch::BranchTree::initLCOs, HPX_NULL, 0);
    }, 0, neuronsCount, NULL);

    compareAllBranches(neuronsCount);

    message("neurox::Branch::finitialize...\n");
    hpx_par_for_sync( [&] (int i, void*) -> int
    {  return hpx_call_sync(neurox::neurons->at(i), Branch::finitialize, HPX_NULL, 0);
    }, 0, neuronsCount, NULL);

#if !defined(NDEBUG) && defined(CORENEURON_H)
    hpx_bcast_rsync(neurox::Input::Coreneuron::Debugger::finitialize);
#endif
    compareAllBranches(neuronsCount);

    int totalSteps = (inputParams->tstop - inputParams->tstart) / inputParams->dt;
    printf("neurox::Algorithm::%s (%d neurons, t=%.03f secs, dt=%.03f milisecs, %d steps):\n",
            inputParams->algorithm == Algorithm::BackwardEulerFixedCommBarrier    ? "BackwardEulerFixedCommBarrier" :
           (inputParams->algorithm == Algorithm::BackwardEulerSlidingTimeWindow   ? "BackwardEulerSlidingTimeWindow" :
           (inputParams->algorithm == Algorithm::BackwardEulerWithPairwiseSteping ? "BackwardEulerWithPairwiseSteping" : "")),
            neuronsCount, inputParams->tstop/1000, inputParams->dt, totalSteps);
    fflush(stdout);

    if (inputParams->algorithm == Algorithm::BackwardEulerSlidingTimeWindow)
    {
        message("neurox::Algorithm::BackwardEulerSlidingTimeWindow::init...\n");
        //hpx_t allreduce = hpx_process_collective_allreduce_new(sizeof(int), _reset, _op);
    }

    hpx_t mainLCO = hpx_lco_and_new(neuronsCount);
    hpx_time_t now = hpx_time_now();
    if (inputParams->algorithm == Algorithm::BackwardEulerFixedCommBarrier)
    {
      for (double t = inputParams->tstart;
         t < inputParams->tstop - inputParams->dt*0.5;
         t += inputParams->dt*Neuron::CommunicationBarrier::commStepSize)
      {
          #ifdef NEUROX_TIME_STEPPING_VERBOSE
            printf("-- t=%.1f ms\n",t);
            fflush(stdout);
          #endif

          hpx_bcast_rsync(neurox::Branch::backwardEulerGroup,
                          &Neuron::CommunicationBarrier::commStepSize, sizeof(int));

          /*TODO compare with this approach
          for (int i=0; i<neuronsCount; i++)
          hpx_call(neurox::neurons->at(i), Branch::backwardEuler, mainLCO,
                   &Neuron::CommunicationBarrier::commStepSize, sizeof(int));
          hpx_lco_wait_reset(mainLCO);
          */

          #if !defined(NDEBUG) && defined(CORENEURON_H)
            if (inputParams->parallelDataLoading) //if parallel execution... spike exchange
              hpx_bcast_rsync(neurox::Input::Coreneuron::Debugger::nrnSpikeExchange);
          #endif
      }
    }
    else if (inputParams->algorithm == Algorithm::BackwardEulerSlidingTimeWindow)
    {
      for (double t = inputParams->tstart;
         t < inputParams->tstop - inputParams->dt*0.5;
         t += inputParams->dt*Neuron::CommunicationBarrier::commStepSize)
      {
        #ifdef NEUROX_TIME_STEPPING_VERBOSE
          printf("-- t=%.1f ms\n",t);
          fflush(stdout);
        #endif

        assert(0);

        for (int i=0; i<neuronsCount; i++)
          hpx_call(neurox::neurons->at(i), Branch::backwardEuler, mainLCO,
                   &Neuron::CommunicationBarrier::commStepSize, sizeof(int));
        hpx_lco_wait_reset(mainLCO);

        #if !defined(NDEBUG) && defined(CORENEURON_H)
          if (inputParams->parallelDataLoading) //if parallel execution... spike exchange
            hpx_bcast_rsync(neurox::Input::Coreneuron::Debugger::nrnSpikeExchange);
        #endif
      }
    }
    else if (inputParams->algorithm == Algorithm::BackwardEulerWithPairwiseSteping)
    {   //NOTE: only comparable at SMP mode or final comparison of results only (due to spike exchange)
        for (int i=0; i<neuronsCount; i++)
          hpx_call(neurox::neurons->at(i), Branch::backwardEuler, mainLCO, &totalSteps, sizeof(int));
        hpx_lco_wait_reset(mainLCO);
    }

    double elapsed = hpx_time_elapsed_ms(now)/1e3;
    compareAllBranches(neuronsCount);
    printf("neurox::end (%d neurons, biological time: %.3f secs, solver time: %.3f secs).\n",
           neuronsCount, inputParams->tstop/1000.0, elapsed);

#if !defined(NDEBUG) && defined(CORENEURON_H)
    neurox::Input::Coreneuron::DataLoader::cleanData(); //if Coreneuron+debug, data will be cleaned up only now
#endif

    hpx_par_for_sync( [&] (int i, void*) -> int
    {  return hpx_call_sync(neurox::neurons->at(i), Branch::clear, HPX_NULL, 0);
    }, 0, neuronsCount, NULL);
    hpx_bcast_rsync(neurox::clear);
    hpx_exit(0,NULL);
}

hpx_action_t clear = 0;
int clear_handler()
{
    neurox_hpx_pin(uint64_t);
    delete [] mechanisms;
    if (myNeurons) { myNeurons->clear(); delete myNeurons; myNeurons = nullptr; }
    if (neurons)   { neurons->clear();   delete neurons;   neurons = nullptr;   }
    neurox_hpx_unpin;
}

void registerHpxActions()
{
    neurox_hpx_register_action(0,neurox::main);
    neurox_hpx_register_action(0,neurox::clear);
    neurox_hpx_register_action(2,neurox::setMechanisms);
    neurox_hpx_register_action(2,neurox::setMechanismsGlobalVars);
}

}; //namespace neurox
