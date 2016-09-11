#include "neurox/neurox.h"
#include "string.h"
#include <numeric>

using namespace neurox;
using namespace neurox::Misc;

class Statistics::SizeInfo
{
  public:
    SizeInfo():neuronId(0), morphologies(0), mechanisms(0), synapses(0),
               metadata(0), globalVars(0), compartmentsCount(0),
               branchesCount(0), mechsInstancesCount(0){};

    ~SizeInfo(){};

    gid_t neuronId;
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
        outstream = fopen(string("neurons-memory-consumption.csv").c_str(), "wt");

    fprintf(outstream, "gid,compartments,branches,mechs-instances,total-KB,morphologies-KB,mechanisms-KB,synapses-KB,metadata-KB\n");

    for (int i=0; i<neuronsCount; i++)
    {
        SizeInfo neuronSize;
        hpx_call_sync(getNeuronAddr(i), Statistics::getNeuronSize, &neuronSize, sizeof(neuronSize));
        fprintf(outstream, "%d,%llu,%llu,%llu,%.1f,%.2f,%.2f,%.2f,%.2f\n",
                neuronSize.neuronId, neuronSize.compartmentsCount, neuronSize.branchesCount,
                neuronSize.mechsInstancesCount, neuronSize.getTotalSize(), neuronSize.morphologies,
                neuronSize.mechanisms, neuronSize.synapses, neuronSize.metadata);
        simSize += neuronSize;
    }

    printf(": SUM %llu neurons, %llu branches, %llu compartments, %llu mech instances, %.1f MB\n",
           neuronsCount, simSize.branchesCount, simSize.compartmentsCount,
           simSize.mechsInstancesCount, simSize.getTotalSize()/1024);
    printf(": AVG per neuron: %.2f branches, %.2f compartments, %.2f mech instances, %.2f KB\n",
           simSize.branchesCount / (double) neuronsCount,
           simSize.compartmentsCount / (double) neuronsCount,
           simSize.mechsInstancesCount / (double) neuronsCount,
           simSize.getTotalSize() / (double) neuronsCount);
    printf(": SUM morphologies %.2f MB, mechanisms %.2f MB, synapses %.2f MB, metadata %.2f MB;\n",
           simSize.morphologies/1024., simSize.mechanisms/1024.,
           simSize.synapses/1024, simSize.metadata/1024);
    printf(": AVG per neuron: morphologies %.2f KB, mechanisms %.2f KB, synapses %.2f KB, metadata %.2f KB;\n",
           simSize.morphologies/ (double) neuronsCount,
           simSize.mechanisms  / (double) neuronsCount,
           simSize.synapses    / (double) neuronsCount,
           simSize.metadata    / (double) neuronsCount );
    printf(": Global vars: %.2f KB (Global data %.2f KB * %d localities)\n",
           simSize.globalVars, simSize.globalVars/HPX_LOCALITIES, HPX_LOCALITIES);
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
    branchSize.morphologies += local->p ? (double) (local->n* sizeof(offset_t))/1024 : 0;
    if (local->neuronTree)
        branchSize.morphologies += local->neuronTree->branches ? (double) (local->neuronTree->branchesCount*sizeof(hpx_t))/1024 : 0;
    branchSize.metadata += (double) sizeof(Branch)/1024;
    branchSize.metadata += (double) sizeof(Branch::MechanismInstance)*mechanismsCount/1024;

    for (int m=0; m<mechanismsCount; m++)
    {
        if (local->mechsInstances[m].count == 0)
            continue;

        branchSize.mechsInstancesCount += local->mechsInstances[m].count;
        branchSize.mechanisms += (double) (sizeof(offset_t) * local->mechsInstances[m].count) /1024;
        if (mechanisms[m]->dataSize>0)
            branchSize.mechanisms += (double) (sizeof(floble_t) * mechanisms[m]->dataSize * local->mechsInstances[m].count)/1024;
        if (mechanisms[m]->pdataSize>0)
            branchSize.mechanisms += (double) (sizeof(offset_t) * mechanisms[m]->pdataSize * local->mechsInstances[m].count)/1024;

        if (mechanisms[m]->pntMap>0)
        {
            if (mechanisms[m]->type == IClamp)
                branchSize.synapses += (double) (sizeof(void*)*2 + sizeof(Point_process)) /1024;
            if (mechanisms[m]->type == ProbAMPANMDA_EMS || mechanisms[m]->type == ProbGABAAB_EMS)
                branchSize.synapses += (double) (sizeof(void*)*3 + sizeof(Point_process) + sizeof(/*RNG_state*/uint32_t[4])) /1024;
        }
    }
    //call the print function in children branches, pass their size to parent branch
    if (local->neuronTree)
    {
      int branchesCount = local->neuronTree->branchesCount; 
      hpx_t branchesLCO = hpx_lco_and_new(branchesCount);
      SizeInfo * subBranchSizes = new SizeInfo[branchesCount];

      for (int c=0; c<branchesCount; c++)
          hpx_call(local->neuronTree->branches[c], Statistics::getNeuronSize,
                   branchesLCO, &subBranchSizes[c], sizeof(SizeInfo) );

      hpx_lco_wait(branchesLCO);

      for (int c=0; c<branchesCount; c++)
          branchSize+=subBranchSizes[c];

      hpx_lco_delete(branchesLCO, HPX_NULL);
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
    printf(": Total mechs instances: %lld\n", totalMechsInstances);

    FILE *outstream = stdout;
    if (writeToFile)
        outstream = fopen(string("mechs-distribution.csv").c_str(), "wt");
    fprintf(outstream, "mech-type,name,instances,avg-per-neuron\n");

    for (int m=0; m<mechanismsCount; m++)
        fprintf(outstream, "%d,%s,%d,%.2f\n", mechanisms[m]->type, mechanisms[m]->sym,
                sumMechsCountPerType[m],(double)sumMechsCountPerType[m]/neuronsCount);

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
    if (local->neuronTree)
    {
      int branchesCount = local->neuronTree->branchesCount;
      hpx_t branchesLCO = hpx_lco_and_new(branchesCount);
      unsigned ** mechsCountPerTypeChild = new unsigned *[branchesCount];

      for (int c=0; c<branchesCount; c++)
          mechsCountPerTypeChild[c] = new unsigned[mechanismsCount];

      for (int c=0; c<branchesCount; c++)
          hpx_call(local->neuronTree->branches[c], Statistics::getNeuronMechanismsDistribution,
                 branchesLCO, mechsCountPerTypeChild[c], sizeof(unsigned)*mechanismsCount);

      hpx_lco_wait(branchesLCO);

      for (int c=0; c<branchesCount; c++)
          for (int m=0; m<mechanismsCount; m++)
              mechsCountPerType[m] += mechsCountPerTypeChild[c][m];

      hpx_lco_delete(branchesLCO, HPX_NULL);
      for (int c=0; c<branchesCount; c++)
          delete [] mechsCountPerTypeChild[c];
      delete [] mechsCountPerTypeChild;
    }
    neurox_hpx_unpin_continue(mechsCountPerType);
}

void Statistics::registerHpxActions()
{
    neurox_hpx_register_action(0, Statistics::getNeuronSize);
    neurox_hpx_register_action(0, Statistics::getNeuronMechanismsDistribution);
}
