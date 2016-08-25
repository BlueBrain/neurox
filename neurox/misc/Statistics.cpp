#include "neurox/neurox.h"
#include "string.h"
#include <numeric>

using namespace neurox;
using namespace neurox::Misc;

class Statistics::SizeInfo
{
  public:
    SizeInfo():neuronId(-1), morphologies(0), mechanisms(0), synapses(0), metadata(0),
        globalVars(0), compartmentsCount(0), branchesCount(0), mechsInstancesCount(0){};
    ~SizeInfo(){};

    int neuronId;
    double morphologies;
    double mechanisms;
    double synapses;
    double metadata;
    double globalVars;
    unsigned long long compartmentsCount;
    unsigned long long branchesCount;
    unsigned long long mechsInstancesCount;

    double getTotalSize(){return morphologies+mechanisms+synapses+metadata+globalVars;}

    SizeInfo& operator+=(const SizeInfo& rhs) {
      mechanisms   += rhs.mechanisms;
      metadata     += rhs.metadata;
      morphologies += rhs.morphologies;
      synapses     += rhs.synapses;
      globalVars   += rhs.globalVars;
      compartmentsCount += rhs.compartmentsCount;
      branchesCount += rhs.branchesCount;
      mechsInstancesCount += rhs.mechsInstancesCount;
      return *this;
    }
};

void Statistics::printSimulationSize(bool writeToFile)
{
    SizeInfo simSize;
    simSize.globalVars = (double) (sizeof(hpx_t) + sizeof(int)*2 + sizeof(Mechanism)*mechanismsCount
                          + sizeof(neurox::Input::InputParams) * HPX_LOCALITIES) /1024;

    FILE *outstream = stdout;
    if (writeToFile)
        outstream = fopen(string("simulation-memory-consumption.txt").c_str(), "wt");

    fprintf(outstream, "Simulation Global vars: %.2f KB (Global data %.2f KB * %d localities)\n",
           simSize.globalVars, simSize.globalVars/HPX_LOCALITIES, HPX_LOCALITIES);

    for (int i=0; i<neuronsCount; i++)
    {
        SizeInfo neuronSize;
        hpx_call_sync(getNeuronAddr(i), Statistics::getNeuronSize, &neuronSize, sizeof(neuronSize));
        fprintf(outstream, "%d: gid %d, %llu compartments, %llu branches, %llu mech instances, %.0f KB:\n",
                i, neuronSize.neuronId, neuronSize.compartmentsCount, neuronSize.branchesCount, neuronSize.mechsInstancesCount, neuronSize.getTotalSize());
        fprintf(outstream, "    morphologies %.2f KB, mechanisms %.2f KB, synapses %.2f KB, metadata %.2f KB\n",
                neuronSize.morphologies, neuronSize.mechanisms, neuronSize.synapses, neuronSize.metadata);
        simSize += neuronSize;
    }
    fprintf(outstream, "TOTAL: %llu neurons, %llu branches, %llu compartments, %llu mech instances, %.0f KB:\n",
           neuronsCount, simSize.branchesCount, simSize.compartmentsCount, simSize.mechsInstancesCount, simSize.getTotalSize());
    fprintf(outstream, "       morphologies %.0f KB, mechanisms %.0f KB, synapses %.0f KB, metadata %.0f KB, globalVars %.0f KB\n",
           simSize.morphologies, simSize.mechanisms, simSize.synapses, simSize.metadata, simSize.globalVars);

    if (writeToFile)
        fclose(outstream);
}

hpx_action_t Statistics::getNeuronSize=0;
int Statistics::getNeuronSize_handler()
{
    neurox_hpx_pin(Branch);
    assert(local->n>0);
    SizeInfo branchSize;
    if (local->soma)
    {
        branchSize.neuronId += local->soma->id;
        branchSize.metadata += (double) sizeof(Neuron) / 1024;
        branchSize.synapses += (double) (local->soma->getNetConsCount()*sizeof(hpx_t)) /1024;
    }
    branchSize.branchesCount++;
    branchSize.compartmentsCount+=local->n;
    branchSize.morphologies += (double) (local->n*(sizeof(floble_t)*6))/1024; //a,b,d,v,rhs,area
    branchSize.morphologies += local->p ? (double) (local->n* sizeof(index_t))/1024 : 0;
    branchSize.morphologies += local->neuronTree->branches ? (double) (local->neuronTree->branchesCount*sizeof(hpx_t))/1024 : 0;
    branchSize.metadata += (double) sizeof(Branch)/1024;
    branchSize.metadata += (double) sizeof(Branch::MechanismInstance)*mechanismsCount/1024;

    for (int m=0; m<mechanismsCount; m++)
    {
        if (local->mechsInstances[m].count == 0)
            continue;

        branchSize.mechsInstancesCount += local->mechsInstances[m].count;
        branchSize.mechanisms += (double) (sizeof(index_t) * local->mechsInstances[m].count) /1024;
        if (mechanisms[m]->dataSize>0)
            branchSize.mechanisms += (double) (sizeof(floble_t) * mechanisms[m]->dataSize * local->mechsInstances[m].count)/1024;
        if (mechanisms[m]->pdataSize>0)
            branchSize.mechanisms += (double) (sizeof(index_t) * mechanisms[m]->pdataSize * local->mechsInstances[m].count)/1024;

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
    for (int c=0; c<local->neuronTree->branchesCount; c++)
    {
        SizeInfo subBranchSize;
        hpx_call_sync(local->neuronTree->branches[c], Statistics::getNeuronSize, &subBranchSize, sizeof(subBranchSize));
        branchSize+=subBranchSize;
    }
    neurox_hpx_unpin_continue(branchSize);
}

void Statistics::printMechanismsDistribution(bool writeToFile)
{
    vector<unsigned> sumMechsCountPerType(mechanismsCount);
    unsigned long long totalMechsInstances=0;
    for (int i=0; i<neuronsCount; i++)
    {
        unsigned mechsCountPerType[mechanismsCount];
        hpx_call_sync(getNeuronAddr(i), Statistics::getNeuronMechanismsDistribution,
                      mechsCountPerType, sizeof(mechsCountPerType));
        for (int m=0; m<mechanismsCount; m++)
        {
            sumMechsCountPerType[m] += mechsCountPerType[m];
            totalMechsInstances += mechsCountPerType[m];
        }
    }

    FILE *outstream = stdout;
    if (writeToFile)
        outstream = fopen(string("mechs-distribution.txt").c_str(), "wt");

    fprintf(outstream, "Total mechs instances: %lld\n", totalMechsInstances);
    for (int m=0; m<mechanismsCount; m++)
        fprintf(outstream, "- type %03d: %d instances of %s (%d), avg per neuron %.1f\n",
             mechanisms[m]->type, sumMechsCountPerType[m], mechanisms[m]->sym, mechanisms[m]->type,
             sumMechsCountPerType[m], (double) sumMechsCountPerType[m]/neuronsCount);

    if (writeToFile)
        fclose(outstream);
}

hpx_action_t Statistics::getNeuronMechanismsDistribution=0;
int Statistics::getNeuronMechanismsDistribution_handler()
{
    neurox_hpx_pin(Branch);
    unsigned mechsCountPerType[mechanismsCount];
    for (int m=0; m<mechanismsCount; m++)
        mechsCountPerType[m] = local->mechsInstances[m].count;

    //call the print function in children branches, pass their size to parent branch
    for (int c=0; c<local->neuronTree->branchesCount; c++)
    {
        vector<unsigned> mechsCountPerTypeChild (mechanismsCount);
        hpx_call_sync(local->neuronTree->branches[c], Statistics::getNeuronMechanismsDistribution,
                  mechsCountPerTypeChild.data(), sizeof(mechsCountPerType));
        for (int m=0; m<mechanismsCount; m++)
            mechsCountPerType[m] += mechsCountPerTypeChild[m];
    }
    neurox_hpx_unpin_continue(mechsCountPerType);
}

void Statistics::registerHpxActions()
{
    neurox_hpx_register_action(0, Statistics::getNeuronSize);
    neurox_hpx_register_action(0, Statistics::getNeuronMechanismsDistribution);
}
