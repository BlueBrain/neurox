#include "neurox/Neurox.h"
#include "string.h"

using namespace NeuroX;
using namespace NeuroX::Misc;

class Statistics::SizeInfo
{
  public:
    SizeInfo():morphologies(0), mechanisms(0), synapses(0), metadata(0), globalVars(0){};
    ~SizeInfo(){};

    double morphologies;
    double mechanisms;
    double synapses;
    double metadata;
    double globalVars;

    double getTotal(){return morphologies+mechanisms+synapses+metadata+globalVars;}

    SizeInfo& operator+=(const SizeInfo& rhs) {
      mechanisms   += rhs.mechanisms;
      metadata     += rhs.metadata;
      morphologies += rhs.morphologies;
      synapses     += rhs.synapses;
      globalVars   += rhs.globalVars;
      return *this;
    }
};

void Statistics::printSimulationSize()
{
    SizeInfo simSize;
    simSize.globalVars = (double) (sizeof(hpx_t) + sizeof(int)*2 + sizeof(Mechanism)*mechanismsCount
                          + sizeof(NeuroX::Input::InputParams) * HPX_LOCALITIES) /1024;
    //printf("Simulation Global vars: %.2f KB (Global data %.2f KB * %d localities)\n",
    //       simSize.globalVars, simSize.globalVars/HPX_LOCALITIES, HPX_LOCALITIES);

    for (int i=0; i<neuronsCount; i++)
    {
        SizeInfo neuronSize;
        hpx_call_sync(getNeuronAddr(i), Statistics::printNeuronSize, &neuronSize, sizeof(neuronSize));
        fflush(stdout);
        simSize += neuronSize;
    }
    printf("Simulation: %.2f KB (morphologies %.2f KB, mechanisms %.2f KB, synapses %.2f KB, metadata %.2f KB, globalVars %.2f KB)\n",
           simSize.getTotal(), simSize.morphologies, simSize.mechanisms, simSize.synapses, simSize.metadata, simSize.globalVars);
}

hpx_action_t Statistics::printNeuronSize=0;
int Statistics::printNeuronSize_handler()
{
    neurox_hpx_pin(Neuron);
    SizeInfo neuronSize;
    neuronSize.metadata = (double) sizeof(Neuron) / 1024;
    neuronSize.synapses = (double) (local->synapses.size()*sizeof(hpx_t)) /1024;
    SizeInfo somaSize;
    hpx_call_sync(local->soma, Statistics::printBranchSize, &somaSize, sizeof(somaSize));
    neuronSize += somaSize;
    printf("  - Neuron %lld: %.0f KB (morphologies %.0f KB, mechanisms %.2\0f KB, synapses %.0f KB, metadata %.0f KB)\n",
           local->id, neuronSize.getTotal(), neuronSize.morphologies, neuronSize.mechanisms, neuronSize.synapses, neuronSize.metadata);
    neurox_hpx_unpin_continue(neuronSize);
}

hpx_action_t Statistics::printBranchSize=0;
int Statistics::printBranchSize_handler()
{
    neurox_hpx_pin(Branch);
    assert(local->n>0);
    SizeInfo branchSize;
    branchSize.morphologies += (double) (local->n*(sizeof(floble_t)*6))/1024; //a,b,d,v,rhs,area
    branchSize.morphologies += local->p ? (double) (local->n* sizeof(index_t))/1024 : 0;
    branchSize.morphologies += local->branches ? (double) (local->branchesCount*sizeof(hpx_t))/1024 : 0;
    branchSize.metadata += (double) sizeof(Branch)/1024;
    branchSize.metadata += (double) sizeof(Branch::MechanismInstance)*mechanismsCount/1024;
    for (int m=0; m<mechanismsCount; m++)
    {
        if (local->mechsInstances[m].instancesCount == 0)
            continue;

        branchSize.mechanisms += (double) (sizeof(index_t) * local->mechsInstances[m].instancesCount) /1024;
        if (mechanisms[m]->dataSize>0)
            branchSize.mechanisms += (double) (sizeof(floble_t) * mechanisms[m]->dataSize * local->mechsInstances[m].instancesCount)/1024;
        if (mechanisms[m]->pdataSize>0)
            branchSize.mechanisms += (double) (sizeof(index_t) * mechanisms[m]->pdataSize * local->mechsInstances[m].instancesCount)/1024;

        if (mechanisms[m]->pntMap>0)
        {
            if (mechanisms[m]->type == IClamp)
                branchSize.synapses += (double) (sizeof(void*)*2 + sizeof(Point_process)) /1024;
            if (mechanisms[m]->type == ProbAMPANMDA_EMS || mechanisms[m]->type == ProbGABAAB_EMS)
                branchSize.synapses += (double) (sizeof(void*)*3 + sizeof(Point_process) + sizeof(/*RNG_state*/uint32_t[4])) /1024;
        }
    }
   // printf("  - Branch (%d compartments): total %.2f KB (morphology %.1f KB, mechanisms %.1f KB, metadata %.3f KB)\n",
   //         local->n, branchSize.getTotal(), branchSize.morphologies, branchSize.mechanisms, branchSize.metadata);

    //call the print function in children branches, pass their size to parent branch
    for (int c=0; c<local->branchesCount; c++)
    {
        SizeInfo subBranchSize;
        hpx_call_sync(local->branches[c], Statistics::printBranchSize, &subBranchSize, sizeof(subBranchSize));
        fflush(stdout);
        branchSize+=subBranchSize;
    }
    neurox_hpx_unpin_continue(branchSize);
}

void Statistics::registerHpxActions()
{
    neurox_hpx_register_action(0, Statistics::printNeuronSize);
    neurox_hpx_register_action(0, Statistics::printBranchSize);
}
