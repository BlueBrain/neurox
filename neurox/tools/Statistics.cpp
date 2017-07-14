#include "neurox/neurox.h"
#include <numeric>
#include <string>
#include <sstream>
#include <iostream>
#include "coreneuron/utils/randoms/nrnran123.h" //RNG data structures

using namespace neurox;
using namespace neurox::tools;

class Statistics::SizeInfo
{
  public:
    SizeInfo():neuronId(0), morphologies(0), mechanisms(0), synapses(0),
               metadata(0), globalVars(0), compartmentsCount(0),
               branchesCount(0), mechsInstancesCount(0){};

    ~SizeInfo(){};

    neuron_id_t neuronId;
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

void Statistics::OutputSimulationSize(bool writeToFile)
{
    SizeInfo simSize;
    simSize.globalVars = (double) (sizeof(hpx_t) + sizeof(int)*2 + sizeof(Mechanism)*mechanismsCount
                          + sizeof(neurox::tools::CmdLineParser) * HPX_LOCALITIES) /1024;

    FILE *outstream = stdout;
    if (writeToFile)
        outstream = fopen(string("neurons-memory-consumption.csv").c_str(), "wt");

    fprintf(outstream, "gid,compartments,branches,mechs-instances,total-KB,morphologies-KB,mechanisms-KB,synapses-KB,metadata-KB\n");

    int neuronsCount = neurons->size();
    for (int i=0; i<neuronsCount; i++)
    {
        SizeInfo neuronSize;
        hpx_call_sync(neurox::neurons->at(i), Statistics::GetNeuronSize, &neuronSize, sizeof(neuronSize));
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

hpx_action_t Statistics::GetNeuronSize=0;
int Statistics::GetNeuronSize_handler()
{
    neurox_hpx_pin(Branch);
    assert(local->nt->end>0);
    SizeInfo branchSize;
    int n = local->nt->end;
    if (local->soma)
    {
        branchSize.neuronId += local->soma->gid;
        branchSize.metadata += (double) sizeof(Neuron) / 1024;
        branchSize.synapses += (double) (local->soma->GetSynapsesCount()*sizeof(Neuron::Synapse)) /1024;
    }
    branchSize.branchesCount++;
    branchSize.compartmentsCount += n;
    branchSize.morphologies += (double) (n*(sizeof(floble_t)*6))/1024; //a,b,d,v,rhs,area
    branchSize.morphologies += local->nt->_v_parent_index ? (double) (n* sizeof(offset_t))/1024 : 0;
    if (local->branchTree)
        branchSize.morphologies += local->branchTree->branches ? (double) (local->branchTree->branchesCount*sizeof(hpx_t))/1024 : 0;
    branchSize.metadata += (double) sizeof(Branch)/1024;
    branchSize.metadata += (double) sizeof(Memb_list)*mechanismsCount/1024;

    for (int m=0; m<mechanismsCount; m++)
    {
        if (local->mechsInstances[m].nodecount == 0)
            continue;

        branchSize.mechsInstancesCount += local->mechsInstances[m].nodecount;
        branchSize.mechanisms += (double) (sizeof(offset_t) * local->mechsInstances[m].nodecount) /1024;
        if (mechanisms[m]->dataSize>0)
            branchSize.mechanisms += (double) (sizeof(floble_t) * mechanisms[m]->dataSize * local->mechsInstances[m].nodecount)/1024;
        if (mechanisms[m]->pdataSize>0)
            branchSize.mechanisms += (double) (sizeof(offset_t) * mechanisms[m]->pdataSize * local->mechsInstances[m].nodecount)/1024;
        if (mechanisms[m]->vdataSize>0)
            branchSize.mechanisms += (double) (sizeof(void*)    * mechanisms[m]->vdataSize * local->mechsInstances[m].nodecount)/1024;

        int type = mechanisms[m]->type;
        if (type == IClamp || type == ProbAMPANMDA_EMS || type == ProbGABAAB_EMS)
            branchSize.synapses += ((double) sizeof(Point_process)  ) /1024;
        if (type == StochKv || type == ProbAMPANMDA_EMS || type == ProbGABAAB_EMS)
            branchSize.synapses += ((double) sizeof(nrnran123_State)) /1024;
    }

    //call the print function in children branches, pass their size to parent branch
    if (local->branchTree && local->branchTree->branchesCount>0)
    {
      int branchesCount = local->branchTree->branchesCount;
      SizeInfo * subBranchSizes = new SizeInfo[branchesCount];

      hpx_t * futures = new hpx_t[branchesCount];
      void ** addrs   = new void*[branchesCount];
      size_t* sizes   = new size_t[branchesCount];
      for (offset_t c = 0; c < branchesCount; c++)
      {
          futures[c] = hpx_lco_future_new(sizeof(SizeInfo));
          addrs[c]   = &subBranchSizes[c];
          sizes[c]   = sizeof(SizeInfo);
          hpx_call(local->branchTree->branches[c],
                   Statistics::GetNeuronSize, futures[c]);
      }
      hpx_lco_get_all(branchesCount, futures, sizes, addrs, NULL);
      hpx_lco_delete_all(branchesCount, futures, NULL);

      delete [] futures;
      delete [] addrs;
      delete [] sizes;

      for (int c=0; c<branchesCount; c++)
          branchSize+=subBranchSizes[c];

      delete [] subBranchSizes;
    }
    neurox_hpx_unpin_continue(branchSize);
}

void Statistics::OutputMechanismsDistribution(bool writeToFile)
{
    unsigned * mechsCountPerType = new unsigned[mechanismsCount];
    unsigned * sumMechsCountPerType = new unsigned[mechanismsCount]();
    unsigned long long totalMechsInstances=0;
    for (int i=0; i<neurox::neurons->size(); i++)
    {
        hpx_call_sync(neurox::neurons->at(i), Statistics::GetNeuronMechanismsDistribution,
                      mechsCountPerType, sizeof(unsigned)*mechanismsCount);
        for (int m=0; m<mechanismsCount; m++)
        {
            sumMechsCountPerType[m] += mechsCountPerType[m];
            totalMechsInstances += mechsCountPerType[m];
        }
    }
    delete [] mechsCountPerType;
    printf(": Total mechs instances: %lld\n", totalMechsInstances);

    FILE *outstream = stdout;
    if (writeToFile)
        outstream = fopen(string("mechs-distribution.csv").c_str(), "wt");
    fprintf(outstream, "mech-type,name,instances,avg-per-neuron\n");

    for (int m=0; m<mechanismsCount; m++)
        fprintf(outstream, "%d,%s,%d,%.2f\n", mechanisms[m]->type, mechanisms[m]->sym,
                sumMechsCountPerType[m],(double)sumMechsCountPerType[m]/neurox::neurons->size());

    if (writeToFile)
        fclose(outstream);

    delete [] sumMechsCountPerType;
}

hpx_action_t Statistics::GetNeuronMechanismsDistribution=0;
int Statistics::GetNeuronMechanismsDistribution_handler()
{
    neurox_hpx_pin(Branch);
    unsigned mechsCountPerType[mechanismsCount];
    for (int m=0; m<mechanismsCount; m++)
        mechsCountPerType[m] = local->mechsInstances[m].nodecount;

    //call the function on children branches, pass their size to parent branch
    if (local->branchTree && local->branchTree->branchesCount>0)
    {
      int branchesCount = local->branchTree->branchesCount;

      unsigned ** mechsCountPerTypeChild = new unsigned *[branchesCount];
      for (int c=0; c<branchesCount; c++)
          mechsCountPerTypeChild[c] = new unsigned[mechanismsCount];

      hpx_t * futures = new hpx_t[branchesCount];
      void ** addrs   = new void*[branchesCount];
      size_t* sizes   = new size_t[branchesCount];
      for (offset_t c = 0; c < branchesCount; c++)
      {
          futures[c] = hpx_lco_future_new(sizeof(unsigned[mechanismsCount]));
          addrs[c]   = mechsCountPerTypeChild[c];
          sizes[c]   = sizeof(unsigned[mechanismsCount]);
          hpx_call(local->branchTree->branches[c],
                   Statistics::GetNeuronMechanismsDistribution, futures[c]);
      }
      hpx_lco_get_all(branchesCount, futures, sizes, addrs, NULL);
      hpx_lco_delete_all(branchesCount, futures, NULL);

      delete [] futures;
      delete [] addrs;
      delete [] sizes;

      for (int c=0; c<branchesCount; c++)
          for (int m=0; m<mechanismsCount; m++)
              mechsCountPerType[m] += mechsCountPerTypeChild[c][m];

      for (int c=0; c<branchesCount; c++)
          delete [] mechsCountPerTypeChild[c];
      delete [] mechsCountPerTypeChild;
    }
    neurox_hpx_unpin_continue(mechsCountPerType);
}

void Statistics::RegisterHpxActions()
{
    neurox_hpx_register_action(neurox_zero_var_action, Statistics::GetNeuronSize);
    neurox_hpx_register_action(neurox_zero_var_action, Statistics::GetNeuronMechanismsDistribution);
}
