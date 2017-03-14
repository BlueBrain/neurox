#include "neurox/neurox.h"
#include <cstring>
#include <map>

#include "nrniv/nrniv_decl.h"
#include "nrniv/nrn_stats.h"

namespace neurox
{
int neuronsCount=-1;
hpx_t neuronsAddr = HPX_NULL;

int mechanismsCount=-1;
extern int * mechanismsMap = nullptr;
neurox::Mechanism ** mechanisms = nullptr;
Input::InputParams * inputParams = nullptr;

hpx_t getNeuronAddr(int i) {
    return hpx_addr_add(neuronsAddr, sizeof(Branch)*i, sizeof(Branch));
}

Mechanism * getMechanismFromType(int type) {
    assert(mechanismsMap[type]!=-1);
    return mechanisms[mechanismsMap[type]];
}

hpx_action_t setNeurons = 0;
int setNeurons_handler(const int nargs, const void *args[], const size_t[])
{
    /** nargs=2 where
     * args[0] = neuronsCount
     * args[1] = neuronsAddr
     */

    neurox_hpx_pin(uint64_t);
    assert(nargs==2);
    neuronsCount = *(int*)args[0];
    neuronsAddr = *(hpx_t*)args[1];
    neurox_hpx_unpin;
}

hpx_action_t setInputParams = 0;
int setInputParams_handler(const Input::InputParams * ip, const size_t)
{
    neurox_hpx_pin(uint64_t);
    if (inputParams!=nullptr)
        delete [] neurox::inputParams;

    neurox::inputParams = new neurox::Input::InputParams();
    memcpy(neurox::inputParams, ip, sizeof(neurox::Input::InputParams));
    neurox_hpx_unpin;
}

hpx_action_t setMechanisms = 0;
int setMechanisms_handler(const int nargs, const void *args[], const size_t sizes[])
{
    /**
     * nargs=3 where:
     * args[0] = array of all mechanisms info
     * args[1] = array of all mechanisms dependencies (parents in mechanisms tree)
     * args[2] = array of all mechanisms successors (children in mechanisms tree)
     * args[3] = array of all mechanisms names (sym)
     */

    neurox_hpx_pin(uint64_t);
    assert(nargs==4);
    mechanismsCount = sizes[0]/sizeof(Mechanism);
    mechanisms = new Mechanism*[mechanismsCount];

    int offsetSuccessors=0, offsetDependencies=0;
    int offsetSym=0;
    int maxMechType=-1;
    for (int m=0; m<mechanismsCount; m++)
    {
        Mechanism & mech = ((Mechanism*) args[0])[m];
        int * dependenciesIds = mech.dependenciesCount == 0 ? nullptr : &((int*) args[1])[offsetDependencies];
        int * successorsIds = mech.successorsCount == 0 ? nullptr : &((int*) args[2])[offsetSuccessors];
        char * sym = mech.symLength == 0 ? nullptr : &((char*) args[3])[offsetSym];
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
    mechanismsMap = new int[maxMechType+1];
    for (int i=0; i<maxMechType+1; i++)
        mechanismsMap[i]=-1;
    for (int m=0; m<mechanismsCount; m++)
        mechanismsMap[mechanisms[m]->type]=m;

    //initializes parent ion index
    for (int m=0; m<mechanismsCount; m++)
    {
      Mechanism * mech = mechanisms[m];
      if (inputParams->multiMex)
      {
        for (int d=0; d<mech->dependenciesCount; d++)
        {
          Mechanism * parent = getMechanismFromType(mech->dependencies[d]);
          if (strcmp("SK_E2", mech->sym)==0 && strcmp("ca_ion", parent->sym)==0) continue; //TODO hard coded exception
          if (parent->getIonIndex() < Branch::MechanismsGraph::IonIndex::size_writeable_ions)
              mech->dependencyIonIndex = parent->getIonIndex();
        }
      }
      else //not applicable
      {
          mech->dependencyIonIndex = Branch::MechanismsGraph::IonIndex::no_ion;
      }
    }
    neurox_hpx_unpin;
}

hpx_action_t main = 0;
static int main_handler(char ***argv, size_t argc)
{
    printf("neurox started (localities: %d, threads/locality: %d)\n", HPX_LOCALITIES, HPX_THREADS);

    //parse command line arguments and broadcasts it to other localities
    Input::InputParams * inputParams_local = new Input::InputParams(argc, *argv);
    hpx_bcast_rsync(neurox::setInputParams, inputParams_local, sizeof (Input::InputParams));
    delete inputParams_local; inputParams_local = nullptr;
    assert(neurox::inputParams != nullptr);

    // many steps with large dt so that cells start at their resting potential
    assert(neurox::inputParams->forwardSkip <= 0); //not supported yet

    //reads morphology data
    printf("neurox::Input::Coreneuron::DataLoader::loadData...\n");
    neurox::Input::Coreneuron::DataLoader::loadData(argc, *argv);

    if (neurox::inputParams->outputStatistics)
    {
      printf("neurox::Misc::Statistics::printMechanismsDistribution...\n", neuronsCount);
      Misc::Statistics::printMechanismsDistribution();
      printf("neurox::Misc::Statistics::printSimulationSize...\n", neuronsCount);
      Misc::Statistics::printSimulationSize();
      //hpx_exit(0,NULL);
    }

    printf("neurox::Branch::NeuronTree::initLCOs...\n");
    hpx_par_for_sync( [&] (int i, void*) -> int
    {  return hpx_call_sync(getNeuronAddr(i), Branch::BranchTree::initLCOs, HPX_NULL, 0);
    }, 0, neuronsCount, NULL);

#if !defined(NDEBUG) && defined(CORENEURON_H)
    printf("NDEBUG::Input::CoreNeuron::DataComparison::compareBranch...\n");
    hpx_par_for_sync( [&] (int i, void*) -> int
    {  return hpx_call_sync(getNeuronAddr(i), Input::Coreneuron::Debugger::compareBranch, HPX_NULL, 0);
    }, 0, neuronsCount, NULL);
#endif

    printf("neurox::BackwardEuler::finitialize...\n");
    hpx_par_for_sync( [&] (int i, void*) -> int
    {  return hpx_call_sync(getNeuronAddr(i), Branch::finitialize, HPX_NULL, 0);
    }, 0, neuronsCount, NULL);

#if !defined(NDEBUG) && defined(CORENEURON_H)
    Input::Coreneuron::Debugger::coreNeuronFinitialize();
    printf("NDEBUG::Input::CoreNeuron::DataComparison::compareBranch...\n");
    hpx_par_for_sync( [&] (int i, void*) -> int
    {  return hpx_call_sync(getNeuronAddr(i), Input::Coreneuron::Debugger::compareBranch, HPX_NULL, 0);
    }, 0, neuronsCount, NULL);
#endif

    int totalSteps = (inputParams->tstop - inputParams->tstart) / inputParams->dt;
    printf("neurox::Algorithm::%s (%d neurons, t=%.03f secs, dt=%.03f milisecs, %d steps):\n",
            inputParams->algorithm == Algorithm::BackwardEulerSyncFixedCommStepDebug ? "BackwardEulerSyncFixedCommStepDebug" :
           (inputParams->algorithm == Algorithm::BackwardEulerAsyncFixedCommStep     ? "BackwardEulerAsyncFixedCommStep" :
           (inputParams->algorithm == Algorithm::BackwardEulerWithPairwiseSteping    ? "BackwardEulerWithPairwiseSteping" : "")),
           neurox::neuronsCount, inputParams->tstop/1000, inputParams->dt, totalSteps);

    hpx_t mainLCO = hpx_lco_and_new(neuronsCount);
    hpx_time_t now = hpx_time_now();
    if (inputParams->algorithm == Algorithm::BackwardEulerSyncFixedCommStepDebug
     || inputParams->algorithm == Algorithm::BackwardEulerAsyncFixedCommStep)
    {
      for (double t = inputParams->tstart;
         t < inputParams->tstop - inputParams->dt*0.5;
         t += inputParams->dt*Neuron::commStepSize)
      {
        printf("-- t=%.1f ms\n",t);
        for (int i=0; i<neuronsCount; i++)
          hpx_call(getNeuronAddr(i), Branch::backwardEuler, mainLCO, &Neuron::commStepSize, sizeof(int));
        hpx_lco_wait_reset(mainLCO);

        //Uncomment to do all neurons comparison at every comm time step
        /*#if !defined(NDEBUG) && defined(CORENEURON_H)
          for (int s=0; s<commStepSize; s++)
            Input::Coreneuron::Debugger::fixed_step_minimal();
          hpx_par_for_sync( [&] (int i, void*) -> int
          {  return hpx_call_sync(getNeuronAddr(i), Input::Coreneuron::Debugger::compareBranch, HPX_NULL, 0);
          }, 0, neuronsCount, NULL);
        #endif */
      }
    }
    else if (inputParams->algorithm == Algorithm::BackwardEulerWithPairwiseSteping)
    {
        for (int i=0; i<neuronsCount; i++)
          hpx_call(getNeuronAddr(i), Branch::backwardEuler, mainLCO, &totalSteps, sizeof(int));
        hpx_lco_wait_reset(mainLCO);

        //Uncomment to compare the final result
        /*#if !defined(NDEBUG) && defined(CORENEURON_H)
          printf("NDEBUG::Input::CoreNeuron::DataComparison::compareBranch...\n");
          for (int s=0; s<totalSteps; s++)
            Input::Coreneuron::Debugger::fixed_step_minimal();
            hpx_par_for_sync( [&] (int i, void*) -> int
              {  return hpx_call_sync(getNeuronAddr(i), Input::Coreneuron::Debugger::compareBranch, HPX_NULL, 0);
              }, 0, neuronsCount, NULL);
        #endif */
    }

    double elapsed = hpx_time_elapsed_ms(now)/1e3;
    printf("neurox::end (%d neurons, biological time: %.3f secs, solver time: %.3f secs).\n",
           neurox::neuronsCount, inputParams->tstop/1000.0, elapsed);

#if !defined(NDEBUG) && defined(CORENEURON_H)
    neurox::Input::Coreneuron::DataLoader::cleanData(); //if Coreneuron+debug, data will be cleaned up only now
#endif
    hpx_exit(0,NULL);
}

hpx_action_t clear = 0;
int clear_handler()
{
    delete [] mechanisms;
}

void registerHpxActions()
{
    neurox_hpx_register_action(1,neurox::main);
    neurox_hpx_register_action(0,neurox::clear);
    neurox_hpx_register_action(2,neurox::setNeurons);
    neurox_hpx_register_action(1,neurox::setInputParams);
    neurox_hpx_register_action(2,neurox::setMechanisms);
}

}; //namespace neurox
