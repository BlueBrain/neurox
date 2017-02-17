#include <stdio.h>
#include <string>
#include <set>
#include <map>
#include <list>
#include <tuple>
#include <algorithm>
#include <numeric>

#include "coreneuron/nrnconf.h"
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrnoc/nrnoc_decl.h"
#include "coreneuron/nrniv/netcon.h"
#include "coreneuron/nrniv/netcvode.h"
#include "coreneuron/nrniv/nrniv_decl.h"
#include "coreneuron/nrniv/nrn_assert.h"
#include "coreneuron/nrniv/nrn_setup.h"
#include "coreneuron/utils/memory_utils.h"
#include "coreneuron/nrnoc/nrnoc_decl.h" //nrn_is_ion()
#include "coreneuron/nrniv/vrecitem.h"

#include "neurox/neurox.h"

using namespace std;
using namespace neurox::Input;
using namespace neurox::Input::Coreneuron;

FILE *fileNetcons;
std::set<neuron_id_t> allNeuronsIdsSet;

neuron_id_t DataLoader::getNeuronIdFromNrnThreadId(int nrn_id)
{
    return (neuron_id_t) nrn_threads[nrn_id].presyns[0].gid_;
}

void DataLoader::printSubClustersToFile(FILE * fileCompartments, Compartment *topCompartment)
{
  if (inputParams->outputCompartmentsDot)
  {
    assert(topCompartment!=NULL);
    fprintf(fileCompartments, "subgraph cluster_%d {\n", topCompartment->id);
    fprintf(fileCompartments, "style=filled; color=blue; fillcolor=floralwhite; node [style=filled, color=floralwhite];\n");
    Compartment * comp=NULL;
    for (comp = topCompartment; comp->branches.size()==1; comp = comp->branches.front())
        //fprintf(fileCompartments, "%d [color=black; fillcolor=white];\n", comp->id);
        fprintf(fileCompartments, "%d [color=gray; fillcolor=white];\n", comp->id);
    //fprintf(fileCompartments, "%d [color=black; fillcolor=white]\n", comp->id);
    fprintf(fileCompartments, "%d [color=gray; fillcolor=white];\n", comp->id);
    fprintf(fileCompartments, "}\n");
    for (int c=0; c<comp->branches.size(); c++)
        printSubClustersToFile(fileCompartments, comp->branches[c]);
  }
}

PointProcInfo DataLoader::getPointProcInfoFromDataPointer(NrnThread * nt, double *pd)
{
    PointProcInfo ppi;
    ppi.nodeId=-1;
    ppi.size = -1; //will be set outside this function
    bool found=false;
    for (NrnThreadMembList* tml = nt->tml; !found && tml!=NULL; tml = tml->next) //For every mechanism
    {
        int type = tml->index;
        Mechanism * mech = getMechanismFromType(type);
        Memb_list * ml = tml->ml;
        unsigned dataOffset=0;
        for (int n=0; n<ml->nodecount; n++)
        {
           // if is this mechanism and this instance
           if (&ml->data[dataOffset] <= pd && pd < &ml->data[dataOffset+mech->dataSize])
           {
               assert(n < 2^sizeof(offset_t));
               ppi.nodeId = ml->nodeindices[n];
               ppi.mechType = type;
               ppi.mechInstance = (offset_t) n;
               ppi.instanceDataOffset = (offset_t) (pd - &ml->data[dataOffset]);
               found=true;
               break;
           }
           dataOffset += (unsigned) mech->dataSize;
        }
    }
    assert(found);
    return ppi;
}

hpx_action_t DataLoader::createNeuron = 0;
int DataLoader::createNeuron_handler(const int *nrnThreadId_ptr, const size_t)
{
        neurox_hpx_pin(Branch);
        NrnThread & nt = nrn_threads[*nrnThreadId_ptr];
        assert(nt.id == *nrnThreadId_ptr);
        neuron_id_t neuronId = getNeuronIdFromNrnThreadId(nt.id);

        #if NETCONS_DOT_OUTPUT_COMPARTMENTS_WITHOUT_NETCONS==true
        if (inputParams->outputNetconsDot)
            fprintf(fileNetcons, "%d [style=filled, shape=ellipse];\n", neuronId);
        #endif

        //======= 1 - reconstructs matrix (solver) values =======

        deque<Compartment*> compartments;
        for (int n=0; n<nt.end; n++)
            compartments.push_back( new Compartment(
                (offset_t) n, (floble_t) nt._actual_a[n],
                (floble_t) nt._actual_b[n], (floble_t) nt._actual_d[n],
                (floble_t) nt._actual_v[n], (floble_t) nt._actual_rhs[n],
                (floble_t) nt._actual_area[n],(offset_t) nt._v_parent_index[n]));

        //reconstructs parents tree
        for (int n=1; n<nt.end; n++) //exclude top (no parent)
        {
            Compartment * parentCompartment = compartments.at(nt._v_parent_index[n]);
            parentCompartment->addChild(compartments.at(n));
        }

        if (inputParams->outputCompartmentsDot)
        {
          FILE *fileCompartments = fopen(string("compartments_"+to_string(neuronId)+".dot").c_str(), "wt");
          fprintf(fileCompartments, "graph G_%d\n{ bgcolor=%s;  node [shape=cylinder];\n", neuronId, DOT_PNG_BACKGROUND_COLOR );
          printSubClustersToFile(fileCompartments, compartments.at(0)); //add subclusters
          for (auto c : compartments) //draw edges
            for (auto k : c->branches)
                fprintf(fileCompartments, "%d -- %d %s;\n", c->id, k->id, "");
          fprintf(fileCompartments, "}\n");
          fclose(fileCompartments);
        }

        //======= 2 - reconstructs mechanisms instances ========

        unsigned vdataTotalOffset=0;
        unsigned dataTotalOffset=nt.end*6;
        unsigned pointProcTotalOffset=0;
        map<offset_t, pair<int,offset_t>> offsetToInstance; //map of data offset -> mech instance (mech Id, node Id)
        for (NrnThreadMembList* tml = nt.tml; tml!=NULL; tml = tml->next) //For every mechanism
        {
            int pdataOffset = 0;
            int dataOffset  = 0;
            int type = tml->index;
            Memb_list * ml = tml->ml; //Mechanisms application to each compartment
            Mechanism * mech = getMechanismFromType(type);
            for (int n=0; n<ml->nodecount; n++) //for every mech instance (or compartment this mech is applied to)
            {
                assert (ml->nodeindices[n] < compartments.size());
                assert (ml->nodeindices[n] < 2^(sizeof(offset_t)*8));
                if (mech->isIon)
                    offsetToInstance[dataTotalOffset] =
                            make_pair(type, (offset_t) ml->nodeindices[n]);
                double * data = &ml->data[dataOffset];
                int * pdata = &ml->pdata[pdataOffset];
                void ** vdata = &nt._vdata[vdataTotalOffset];
                Compartment * compartment = compartments.at(ml->nodeindices[n]);
                assert(compartment->id == ml->nodeindices[n]);
                for (int i=0; i<mech->dataSize; i++)
                {   assert(data[i]==nt._data[dataTotalOffset+i]); }

                if (mech->pntMap > 0) //vdata
                {
                    assert((type==IClamp && mech->pdataSize==2)
                        ||((type==ProbAMPANMDA_EMS || type==ProbGABAAB_EMS) && mech->pdataSize==3));

                    //Offsets in pdata: 0: data (area), 1: point proc in nt._vdata, [2: rng in nt._vdata]
                    for (int v=0; v<mech->vdataSize; v++)
                        compartment->vdata.push_back(vdata[v]);
#ifndef NDEBUG
                    Point_process * pp = &nt.pntprocs[pointProcTotalOffset++];
                    assert(nt._vdata[pdata[1]] == pp);
                    if (mech->pdataSize > 2)
                    {
                        assert(pdata[1]+1 ==pdata[2]); //TODO no need to store offsets 1 and 2 if they are sequential
                        void * RNG = nt._vdata[pdata[2]];
                        (void) RNG;
                    }
                }
                else
                {
                    assert(mech->vdataSize==0);
#endif
                }
                compartment->addMechanismInstance(type, data, mech->dataSize, pdata,  mech->pdataSize);
                dataOffset       += (unsigned) mech->dataSize;
                pdataOffset      += (unsigned) mech->pdataSize;
                dataTotalOffset  += (unsigned) mech->dataSize;
                vdataTotalOffset += (unsigned) mech->vdataSize;
                assert(dataTotalOffset  < 2^sizeof(offsetToInstance));
                assert(dataOffset       < 2^sizeof(offsetToInstance));
                assert(pdataOffset      < 2^sizeof(offsetToInstance));
                assert(vdataTotalOffset < 2^sizeof(offsetToInstance));
            }
        }

        //======= 3 - reconstruct NetCons =====================

        map< neuron_id_t, vector<NetConX*> > netcons ; //netcons per pre-synaptic neuron offset (not id)
        for (int n = 0; n < nt.n_netcon; ++n) {
            NetCon* nc = nt.netcons + n;
            assert(netcon_srcgid.size()>0); //if size==0, then setup_cleanup() in nrn_setup.cpp was called
            assert(nt.id < netcon_srcgid.size());
            assert(netcon_srcgid.at(nt.id)!=NULL);
            int srcGid = netcon_srcgid[nt.id][n];
            assert(srcGid>=0); //gids can be negative! This is a reminder that I should double-check when they happen
            int mechType = nc->target_->_type;
            int weightsCount = pnt_receive_size[mechType];
            size_t weightsOffset = nc->u.weight_index_;
            assert(weightsOffset < nt.n_weight);
            floble_t * weights = &(nt.weights[weightsOffset]);
            netcons[srcGid].push_back(
                new NetConX(mechType, (offset_t) nc->target_->_i_instance,
                         (floble_t) nc->delay_, weights, weightsCount, nc->active_));
        }

        if (inputParams->outputNetconsDot)
        {
          int netConsFromOthers=0;
          for (auto nc : netcons)
          {
            int srcGid = nc.first;
            if (allNeuronsIdsSet.find(srcGid) == allNeuronsIdsSet.end())
                netConsFromOthers++;
            else
            {
                floble_t minDelay=99999;
                for (auto ncv : nc.second) //get minimum delay between neurons
                    minDelay = std::min(minDelay, ncv->delay);
                fprintf(fileNetcons, "%d -> %d [label=\"%d (%.2fms)\"];\n",
                        srcGid , neuronId, nc.second.size(), minDelay);
            }
          }
          #if NETCONS_DOT_OUTPUT_NETCONS_FROM_EXTERNAL_NEURONS==true
            if (netConsFromOthers>0)
              fprintf(fileNetcons, "%s -> %d [label=\"%d\" fontcolor=gray color=gray arrowhead=vee fontsize=12];\n", "external", neuronId, netConsFromOthers);
          #endif
        }

        //======= 4 - reconstruct VecPlayContinuous events =======
        for (int v=0; v<nt.n_vecplay; v++)
        {
            VecPlayContinuous *vec = (VecPlayContinuous*) nt._vecplay[v];
            //discover node, mechanism and data offset id that *pd points to
            PointProcInfo ppi = getPointProcInfoFromDataPointer(&nt, vec->pd_);
            ppi.size = vec->y_->size();
            compartments.at(ppi.nodeId)->addVecPlay(vec->t_->data(), vec->y_->data(), ppi);
        }

        //======= 5 - recursively create branches tree ===========

        floble_t APthreshold = (floble_t) nrn_threads[nt.id].presyns[0].threshold_;
        int thvar_index = nrn_threads[nt.id].presyns[0].thvar_index_;
        createBranch(target, compartments, compartments.at(0), netcons,
                     (size_t) compartments.size(), offsetToInstance);
#ifdef CORENEURON_H
        hpx_call_sync(target, Branch::initSoma, NULL, 0,
                      &neuronId, sizeof(neuron_id_t), &APthreshold, sizeof(floble_t),
                      &thvar_index, sizeof(thvar_index), &nt.id, sizeof(nt.id));
#else
        hpx_call_sync(target, Branch::initSoma, NULL, 0,
                      &neuronId, sizeof(neuron_id_t), &APthreshold, sizeof(floble_t));
#endif
        for (auto c : compartments)
            delete c;
        for (auto nc: netcons)
            for (auto nvc : nc.second)
                delete nvc;
    neurox_hpx_unpin;
}

void DataLoader::cleanData()
{
    nrn_cleanup();
}

void DataLoader::loadData(int argc, char ** argv)
{
    assert (neurox::inputParams!=nullptr); //TODO mem corruption

    cn_input_params input_params;
    nrn_init_and_load_data(argc, argv, input_params, false);

    int neuronsCount = std::accumulate(nrn_threads, nrn_threads+nrn_nthread, 0,
                                 [](int n, NrnThread & nt){return n+nt.ncell;});

    allNeuronsIdsSet = std::set<neuron_id_t>();
    for (int i=0; i<neuronsCount; i++)
    {
        assert(nrn_threads[i].ncell == 1);
        allNeuronsIdsSet.insert(getNeuronIdFromNrnThreadId(i));
    }

    #if COMPARTMENTS_DOT_OUTPUT_CORENEURON_STRUCTURE == true
    if (inputParams->outputCompartmentsDot)
    {
      for (int i=0; i<neuronsCount; i++)
      {
        neuron_id_t neuronId = getNeuronIdFromNrnThreadId(i);
        FILE *fileCompartments = fopen(string("compartments_"+to_string(neuronId)+"_NrnThread.dot").c_str(), "wt");
        fprintf(fileCompartments, "graph G%d\n{  node [shape=cylinder];\n", neuronId );

        //for all nodes in this NrnThread
        NrnThread * nt = &nrn_threads[i];
        for (int n=nt->ncell; n<nt->end; n++)
            fprintf(fileCompartments, "%d -- %d;\n", nt->_v_parent_index[n], n);
        fprintf(fileCompartments, "}\n");
        fclose(fileCompartments);
      }
    }
    #endif

    /** Reconstructs unique data related to each mechanism*
     * nargs=3 where:
     * args[0] = array of all mechanisms info
     * args[1] = array of all mechanisms depedencies (children on mechanisms execution tree)
     * args[2] = array of all mechanisms names (sym)
     */
    std::vector<Mechanism> mechsData;
    std::vector<int> mechsSuccessorsId;
    std::vector<char> mechsSym;

    //Different nrn_threads[i] have diff mechanisms sets; we'll get the union of all neurons' mechs
    std::list<NrnThreadMembList*> uniqueMechs; //list of unique mechanisms
    std::set<int> uniqueMechIds; //list of unique mechanism ids

    //insert all mechs of first neuron
    for (NrnThreadMembList* tml = nrn_threads[0].tml; tml!=NULL; tml = tml->next)
    {
        uniqueMechs.push_back(tml);
        uniqueMechIds.insert(tml->index);
    }

    //insert all mechs from other neurons that do not yet exist
    for (int i=1; i<neuronsCount; i++)
        for (NrnThreadMembList* tml = nrn_threads[i].tml; tml->next!=NULL; tml = tml->next)
            if (uniqueMechIds.find(tml->next->index) == uniqueMechIds.end()) //if next mech does not exist
            {   //find correct position in list and insert it there:
                for (auto it = uniqueMechs.begin(); it != uniqueMechs.end(); it++) //...for all existing mechs
                    if ((*it)->index == tml->index)
                    {
                        auto it_next = std::next(it,1);
                        uniqueMechs.insert(it_next,tml->next); //reminder: .insert adds elements in position before iterator
                        //we have it -> it_new -> it_next. Now we will set the value of next pointer
                        //auto it_new = std::prev(it_next,1);
                        //(*it)->next = *it_new; //Do not change next pointers or neuron is incorrect
                        //(*it_new)->next = *it_next;
                        uniqueMechIds.insert(tml->next->index);
                        break;
                    }
            }

    for (auto tml_it = uniqueMechs.begin(); tml_it != uniqueMechs.end(); tml_it++)
    {
        auto & tml = *tml_it;
        int type = tml->index;
        vector<int> successorsIds;
        int dependenciesCount;

        if (inputParams->multiMex)
        {
          std::set<int> dependenciesIds ; //set of 'unique' dependencies (ignores several dependencies between same mechs pair)
          for (auto & tml2 : uniqueMechs)
          {
            int otherType = tml2->index;
            for (int d=0; d<nrn_prop_dparam_size_[otherType]; d++)
            {
              int ptype = memb_func[otherType].dparam_semantics[d];
              if (otherType == type && ptype>0 && ptype<1000)
                dependenciesIds.insert(ptype); //this mech depends on another one
              if (otherType!= type && ptype==type)
                if (std::find(successorsIds.begin(), successorsIds.end(), otherType) == successorsIds.end())
                  successorsIds.push_back(otherType); //other mech depends on this one
            }
          }
          dependenciesCount = dependenciesIds.size();
        }
        else
        {
          //all except second element (the one after capacitance) have 1 dependency
          dependenciesCount = type==CAP ? 0 : 1;
          if (tml->index != uniqueMechs.back()->index)
          {
              auto tml_next_it = std::next(tml_it,1);
              successorsIds.push_back((*tml_next_it)->index);
          }
        }

        int symLength = memb_func[type].sym ? std::strlen(memb_func[type].sym) : 0;
        if (strcmp(memb_func[type].sym, "PatternStim")==0 && inputParams->patternStim[0]=='\0')
                continue; //only load PatternStim if path initialized

        mechsData.push_back(
            Mechanism (type, nrn_prop_param_size_[type], nrn_prop_dparam_size_[type],
                       nrn_is_artificial_[type], pnt_map[type], nrn_is_ion(type),
                       symLength, nullptr, //sym will be serialized below
                       dependenciesCount, successorsIds.size(), nullptr)); //children will be serialized below

        mechsSuccessorsId.insert(mechsSuccessorsId.end(), successorsIds.begin(), successorsIds.end());
        mechsSym.insert(mechsSym.end(), memb_func[type].sym, memb_func[type].sym + symLength);
    }

    if (inputParams->patternStim[0]!='\0') //"initialized"
    {
        assert(0); //not an error: should be initialized already by coreneuron (and above)
        //in the future this should be contidtional (once we get rid of coreneuron data loading)
    }

    printf("neurox::setMechanisms...\n", mechsData.size());
    hpx_bcast_rsync(neurox::setMechanisms,
                    mechsData.data(), sizeof(Mechanism)*mechsData.size(),
                    mechsSuccessorsId.data(), sizeof(int)* mechsSuccessorsId.size(),
                    mechsSym.data(), sizeof(char)*mechsSym.size());

    mechsData.clear(); mechsSuccessorsId.clear(); mechsSym.clear();

#if !defined(NDEBUG) && defined(CORENEURON_H)
    neurox::Input::Coreneuron::Debugger::compareMechanismsFunctionPointers(uniqueMechs);
#endif

    if (inputParams->outputMechanismsDot)
    {
      FILE *fileMechs = fopen(string("mechanisms.dot").c_str(), "wt");
      fprintf(fileMechs, "digraph G\n{ bgcolor=%s;\n", DOT_PNG_BACKGROUND_COLOR);
      fprintf(fileMechs, "graph [ratio=0.3];\n", "start");
      fprintf(fileMechs, "%s [style=filled, shape=Mdiamond, fillcolor=beige];\n", "start");
      fprintf(fileMechs, "%s [style=filled, shape=Mdiamond, fillcolor=beige];\n", "end");
      fprintf(fileMechs, "\"%s (%d)\" [style=filled, fillcolor=beige];\n",
            getMechanismFromType(CAP)->sym, CAP);
      for (int m =0; m< mechanismsCount; m++)
      {
        Mechanism * mech = mechanisms[m];
        if (mech->pntMap > 0) //if is point process make it dotted
            fprintf(fileMechs, "\"%s (%d)\" [style=dotted];\n", mech->sym, mech->type);

        if (mech->dependenciesCount==0 && mech->type!=CAP) //top mechanism
            fprintf(fileMechs, "%s -> \"%s (%d)\";\n", "start", mech->sym, mech->type);

        if (mech->successorsCount==0 && mech->type!= CAP) //bottom mechanism
            fprintf(fileMechs, "\"%s (%d)\" -> %s;\n", mech->sym, mech->type, "end");

        for (int d=0; d<mech->successorsCount; d++)
            fprintf(fileMechs, "\"%s (%d)\" -> \"%s (%d)\";\n",  mech->sym, mech->type,
                    getMechanismFromType(mech->successors[d])->sym, getMechanismFromType(mech->successors[d])->type);
      }
      fprintf(fileMechs, "}\n");
      fclose(fileMechs);
    }

    //allocate HPX memory space for neurons
    printf("neurox::setNeurons...\n");
    hpx_t neuronsAddr = hpx_gas_calloc_cyclic(neuronsCount, sizeof(Branch), NEUROX_HPX_MEM_ALIGNMENT);
    hpx_bcast_rsync(neurox::setNeurons, &neuronsCount, sizeof(int), &neuronsAddr, sizeof(hpx_t));
    assert(neuronsAddr != HPX_NULL);

    if (inputParams->outputNetconsDot)
    {
      assert(HPX_LOCALITIES == 1);
      fileNetcons = fopen(string("netcons.dot").c_str(), "wt");
      //fprintf(fileNetcons, "digraph G\n{ bgcolor=%s; layout=circo;\n", DOT_PNG_BACKGROUND_COLOR);
      fprintf(fileNetcons, "digraph G\n{ bgcolor=%s;\n", DOT_PNG_BACKGROUND_COLOR);
      #if NETCONS_DOT_OUTPUT_NETCONS_FROM_EXTERNAL_NEURONS==true
        fprintf(fileNetcons, "external [color=gray fontcolor=gray];\n");
      #endif
    }

    printf("neurox::createNeurons...\n");
    hpx_par_for_sync( [&] (int i, void*) -> int
    {  return hpx_call_sync(getNeuronAddr(i), DataLoader::createNeuron, NULL, 0, &i, sizeof(i) );
    }, 0, neuronsCount, NULL);

    //all neurons have been created, every branch will inform pre-syn ids that they are connected
    assert(allNeuronsIdsSet.size() == neuronsCount);
    std::vector<neuron_id_t> allNeuronsIdsVec  (neuronsCount);
    std::vector<hpx_t> allNeuronsAddrVec (neuronsCount);
    for (int i=0; i<neuronsCount; i++)
    {
        allNeuronsIdsVec[i]  = getNeuronIdFromNrnThreadId(i);
        allNeuronsAddrVec[i] = getNeuronAddr(i);
    }

    printf("neurox::Neuron::broadcastNetCons...\n", neuronsCount);
    hpx_t lco_neurons = hpx_lco_and_new(neuronsCount);
    for (int i=0; i<neuronsCount; i++)
        hpx_call(getNeuronAddr(i), DataLoader::broadcastNetCons, lco_neurons,
                 allNeuronsIdsVec.data(),  allNeuronsIdsVec.size() *sizeof(neuron_id_t),
                 allNeuronsAddrVec.data(), allNeuronsAddrVec.size()*sizeof(hpx_t));
    hpx_lco_wait(lco_neurons);
    hpx_lco_delete(lco_neurons, NULL);

    if (inputParams->outputNetconsDot)
    {
      fprintf(fileNetcons, "}\n");
      fclose(fileNetcons);
    }

    nrn_setup_cleanup();
}

void DataLoader::getNetConsBranchData(
        deque<Compartment*> & compartments, map<neuron_id_t, vector<NetConX*> > & netcons,
        vector<NetConX> & branchNetCons, vector<neuron_id_t> & branchNetConsPreId,
        vector<floble_t> & branchNetConsArgs)
{
    //get size and allocate it
    int ncCount=0, weightsCount=0;
    for (auto nc : netcons)
        for (auto nc2 : nc.second)
        {
            ncCount++;
            weightsCount += nc2->weightsCount;
        }

    branchNetCons = vector<NetConX> (ncCount);
    branchNetConsPreId = vector<neuron_id_t> (ncCount);
    branchNetConsArgs = vector<floble_t> (ncCount);

    int i=0;
    for (auto nc : netcons)
        for (auto nc2 : nc.second)
        {
           assert(nc2!=NULL);
           assert(nc2->weights != NULL);
           memcpy(&branchNetCons.at(i), nc2, sizeof(NetConX)); //push_back will call destructor!
           memcpy(&branchNetConsPreId.at(i), &nc.first, sizeof(int));
           branchNetConsArgs.insert(branchNetConsArgs.end(), nc2->weights, nc2->weights + nc2->weightsCount);
           branchNetCons[i].weights=NULL; //needed or will be pointing to original netcons data and will call destructor
           i++;
        }
}

void DataLoader::getVecPlayBranchData(
        deque<Compartment*> & compartments, vector<floble_t> & vecPlayTdata,
        vector<floble_t> & vecPlayYdata, vector<PointProcInfo> & vecPlayInfo)
{
    for (auto comp : compartments)
    {
        vecPlayTdata.insert(vecPlayTdata.end(), comp->vecPlayTdata.begin(), comp->vecPlayTdata.end());
        vecPlayYdata.insert(vecPlayYdata.end(), comp->vecPlayYdata.begin(), comp->vecPlayYdata.end());
        vecPlayInfo.insert(vecPlayInfo.end(), comp->vecPlayInfo.begin(), comp->vecPlayInfo.end());
    }
}

int DataLoader::getBranchData(
        deque<Compartment*> & compartments, vector<floble_t> & data,
        vector<offset_t> & pdata, vector<void*> & vdata, vector<offset_t> & p,
        vector<offset_t> & instancesCount, vector<offset_t> & nodesIndices,
        int totalN, map<offset_t, pair<int, offset_t>> & offsetToInstance)
{
    for (auto comp : compartments)
        { assert (comp != NULL); }
    for (auto comp : compartments)
        data.push_back(comp->rhs);
    for (auto comp : compartments)
        data.push_back(comp->d);
    for (auto comp : compartments)
        data.push_back(comp->a);
    for (auto comp : compartments)
        data.push_back(comp->b);
    for (auto comp : compartments)
        data.push_back(comp->v);
    for (auto comp : compartments)
        data.push_back(comp->area);
    for (auto comp: compartments)
        p.push_back(comp->p);

    vector<vector<floble_t>> dataMechs (mechanismsCount);
    vector<vector<offset_t>> pdataMechs (mechanismsCount);
    vector<vector<offset_t>> nodesIndicesMechs (mechanismsCount);
    vector<vector<void*>> vdataMechs (mechanismsCount);

    int n=0;
    vector<offset_t> pdataType; //type of pdata offset per pdata entry
    map< pair<int, offset_t>, offset_t> instanceToOffset; //from pair of < mech type, OLD node id> to offset in NEW representation
    map<offset_t,offset_t> fromNeuronToBranchId; //map of branch id per compartment id
    for (auto comp : compartments)
    {
        int compDataOffset=0;
        int compPdataOffset=0;
        int compVdataOffset=0;
        fromNeuronToBranchId[comp->id] = n;
        for (int m=0; m<comp->mechsTypes.size(); m++) //for all instances
        {
            int mechType = comp->mechsTypes[m];
            int mechOffset = mechanismsMap[mechType];
            assert(mechOffset>=0 && mechOffset<mechanismsCount);
            Mechanism * mech = mechanisms[mechOffset];
            dataMechs[mechOffset].insert(dataMechs[mechOffset].end(), &comp->data[compDataOffset], &comp->data[compDataOffset+mech->dataSize] );
            pdataMechs[mechOffset].insert(pdataMechs[mechOffset].end(), &comp->pdata[compPdataOffset], &comp->pdata[compPdataOffset+mech->pdataSize] );
            vdataMechs[mechOffset].insert(vdataMechs[mechOffset].end(), &comp->vdata[compVdataOffset], &comp->vdata[compVdataOffset+mech->vdataSize] );
            nodesIndicesMechs[mechOffset].push_back(comp->id);
            instancesCount[mechOffset]++;
            compDataOffset  += mech->dataSize;
            compPdataOffset += mech->pdataSize;
            compVdataOffset += mech->vdataSize;
        }
        n++;
    }

    //merge all mechanisms vectors in the final one
    //store the offset of each mechanism data (for later)
    for (int m=0; m<mechanismsCount; m++)
    {
        Mechanism * mech = mechanisms[m];
        for (offset_t n=0; n<instancesCount[m]; n++)
          for (short d=0; d<mech->pdataSize; d++)
            pdataType.push_back(memb_func[mech->type].dparam_semantics[d]);

        int dataOffset=0;
        int pdataOffset=0;
        int vdataOffset=0;
        for (int i=0; i<nodesIndicesMechs[m].size(); i++) //for all instances
        {
            assert(instanceToOffset.find( make_pair(mech->type, nodesIndicesMechs[m][i]) ) == instanceToOffset.end() );
            if (mech->isIon)
                instanceToOffset[ make_pair(mech->type, nodesIndicesMechs[m][i]) ] = data.size();

            data.insert ( data.end(),  &dataMechs[m][dataOffset ],  &dataMechs[m][dataOffset +mech->dataSize ]);
            pdata.insert(pdata.end(), &pdataMechs[m][pdataOffset], &pdataMechs[m][pdataOffset+mech->pdataSize]);
            vdata.insert(vdata.end(), &vdataMechs[m][vdataOffset], &vdataMechs[m][vdataOffset+mech->vdataSize]);
            nodesIndices.push_back(fromNeuronToBranchId[ nodesIndicesMechs[m][i] ]);

            dataOffset  += mech->dataSize;
            pdataOffset += mech->pdataSize;
            vdataOffset += mech->vdataSize;
        }
        dataMechs[m].clear();
        pdataMechs[m].clear();
        vdataMechs[m].clear();
        nodesIndicesMechs[m].clear();
    }

    if (!inputParams->multiSplix)
      //if there are more than one instace of the same ion mech on a node, this fails!
      {assert(instanceToOffset.size() == offsetToInstance.size());}

    //convert all offsets in pdata to the correct ones
    //depends on the pdata offset type (register_mech.c :: hoc_register_dparam_semantics)
    int vdataOffset=0;
    assert(pdataType.size() == pdata.size());
    assert(pdata.size() < 2^(sizeof(offset_t)*8));

    for (size_t i=0; i<pdata.size(); i++)
    {
        offset_t p = pdata.at(i);
        int ptype = pdataType.at(i);
        switch (ptype)
        {
        case -1:  //"area" (6th field)
        {
            assert(p>=totalN*5 && p<totalN<6);
            offset_t oldId = p-totalN*5;
            offset_t newId = fromNeuronToBranchId[oldId];
            pdata[i] = n*5+newId;
            break;
        }
        case -2: //"iontype"
            //do nothing, its a flag (the 'iontype', see nrnoc/eion.c)
            break;
        case -3: //"cvodeieq"
        case -5: //"pointer"
            assert(0); //not used
            break;
        case -4: //"netsend"
        case -6: //"pntproc"
        case -7: //"bbcorepointer"
            pdata[i] = (offset_t) vdataOffset++;
            break;
        default:
            if (ptype>0 && ptype<1000) //name preffixed by '#'
            {   //ptype is the ion (mechanism) type it depends on
                //pdata[i] is an offset in pdata
                std::pair<int, offset_t> mechNodePair;
                offset_t oldOffset = -1;
                for (auto val : offsetToInstance)
                    if (val.first<=pdata[i])
                    {
                        oldOffset = val.first;
                        mechNodePair=val.second;
                    }
                    else break; //already found
                assert(oldOffset!=-1); //make sure was found
                assert(oldOffset>=totalN*6); //also fails if not found
                assert(mechNodePair.first == ptype); //fails if expected type does not agree with the one in the data offset
                assert(pdata[i]-oldOffset < getMechanismFromType(mechNodePair.first)->dataSize); //instance offset is correct
                offset_t newOffset = instanceToOffset.at(mechNodePair);
                assert(newOffset>=n*6);
                pdata[i] = newOffset + (pdata[i]-oldOffset); //'pdata[i]-oldOffset' is the offset on the data vector for that instance
                assert(pdata[i]>=n*6);
            }
            else if (ptype>=1000) //name not preffixed
            {
                // (concentration = ptype-1000;) //do nothing: value of concentration summed with 1000
            }
            else
                throw std::runtime_error("Unknown pdata type %d (FLAG3)\n");
            break;
        }
    }
    return n;
}

hpx_t DataLoader::createBranch(hpx_t target, deque<Compartment*> & compartments, Compartment * topCompartment,
                               map< neuron_id_t, vector<NetConX*> > & netcons, int totalN, map<offset_t, pair<int,offset_t>> & offsetToInstance)
{
    assert(topCompartment!=NULL);
    offset_t n; //number of compartments in branch
    vector<floble_t> data; //compartments info (RHS, D, A, B, V, AREA)*n
    vector<offset_t> pdata; //pointers to data
    vector<offset_t> p; //parent nodes index
    vector<offset_t> instancesCount (mechanismsCount);
    vector<offset_t> nodesIndices;
    vector<hpx_t> branches;
    vector<void*> vdata; //TODO should go away at some point,

    //Vector Play instances
    vector<floble_t> vecPlayT;
    vector<floble_t> vecPlayY;
    vector<PointProcInfo> vecPlayInfo;

    //branch NetCons
    vector<NetConX> branchNetCons;
    vector<neuron_id_t> branchNetConsPreId;
    vector<floble_t> branchNetConsArgs;

    if (inputParams->multiSplix)
    {
        deque<Compartment*> subSection;
        Compartment * comp = NULL;
        for (comp = topCompartment;
             comp->branches.size()==1;
             comp = comp->branches.front())
        {
            subSection.push_back(comp);
        }
        assert(comp!=NULL);
        subSection.push_back(comp); //bottom of a branch (bifurcation or bottom leaf)

        //create sub-section of branch (comp is the bottom compartment of the branch)
        n = getBranchData(subSection, data, pdata, vdata, p, instancesCount, nodesIndices, totalN, offsetToInstance);
        //TODO broken until we manage to translate mechId+instance in neuron to branch!
        //getVecPlayBranchData(compartments, vecPlayT, vecPlayY, vecPlayInfo);
        //getNetConsBranchData(compartments, netcons, branchNetCons, branchNetConsPreId, branchNetConsArgs);

        //recursively create children branches
        for (size_t c=0; c<comp->branches.size(); c++)
            branches.push_back(createBranch(HPX_NULL, compartments, comp->branches[c], netcons, totalN, offsetToInstance));
    }
    else //Flat a la Coreneuron
    {
        n = getBranchData(compartments, data, pdata, vdata, p, instancesCount, nodesIndices, totalN, offsetToInstance);
        getVecPlayBranchData(compartments, vecPlayT, vecPlayY, vecPlayInfo);
        getNetConsBranchData(compartments, netcons, branchNetCons, branchNetConsPreId, branchNetConsArgs);
    }

    //Allocate HPX Branch (top has already been created on main neurons array)
    hpx_t branchAddr = target==HPX_NULL ? hpx_gas_calloc_local(1, sizeof(Branch), NEUROX_HPX_MEM_ALIGNMENT) : target;
    bool multiSplix = inputParams->multiSplix;

    //initiate branch
    hpx_call_sync(branchAddr, Branch::init, NULL, 0,
                  &n, sizeof(offset_t),
                  HPX_NULL, 0, //NOT USED
                  data.size()>0 ? data.data() : nullptr, sizeof(floble_t)*data.size(),
                  pdata.size()>0 ? pdata.data() : nullptr, sizeof(offset_t)*pdata.size(),
                  instancesCount.data(), instancesCount.size()*sizeof(offset_t),
                  nodesIndices.data(), nodesIndices.size()*sizeof(offset_t),
                  multiSplix ? branches.data() : nullptr, multiSplix ? sizeof(hpx_t)*branches.size() : 0,
                  multiSplix ? nullptr : p.data(), multiSplix ? 0 : sizeof(offset_t)*p.size(),
                  vecPlayT.size() > 0 ? vecPlayT.data() : nullptr, sizeof(floble_t)*vecPlayT.size(),
                  vecPlayY.size() > 0 ? vecPlayY.data() : nullptr, sizeof(floble_t)*vecPlayY.size(),
                  vecPlayInfo.size() > 0 ? vecPlayInfo.data() : nullptr, sizeof(PointProcInfo)*vecPlayInfo.size(),
                  branchNetCons.size() > 0 ? branchNetCons.data() : nullptr, sizeof(NetCon)*branchNetCons.size(),
                  branchNetConsPreId.size() > 0 ? branchNetConsPreId.data() : nullptr, sizeof(neuron_id_t)*branchNetConsPreId.size(),
                  branchNetConsArgs.size() > 0 ? branchNetConsArgs.data() : nullptr, sizeof(floble_t)*branchNetConsArgs.size(),
                  vdata.size()>0 ? vdata.data() : nullptr, sizeof(void*)*vdata.size()
                  );

    return branchAddr;
}

hpx_action_t DataLoader::broadcastNetCons = 0;
int DataLoader::broadcastNetCons_handler(const int nargs, const void *args[], const size_t sizes[])
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(DataLoader::broadcastNetCons, args[0], sizes[0], args[1], sizes[1]);
    assert(nargs==2);

    const neuron_id_t * neuronsIds = (const neuron_id_t*) args[0]; //list of all existing neurons ids
    const hpx_t * neuronsAddr = (const hpx_t*) args[1]; //list of all existing neurons addrs (same order as ids)

    //create index of unique pre-neurons ids
    std::set<neuron_id_t> netconsPreNeuronIds;
    for (auto nc : local->netcons)
        netconsPreNeuronIds.insert(nc.first);

    //inform pre-synaptic neurons (once) that we connect
    for (int i=0; i<sizes[0]/sizeof(neuron_id_t); i++) //for all existing neurons
        //if I'm connected to it (ie is not artificial)
        if (netconsPreNeuronIds.find(neuronsIds[i]) != netconsPreNeuronIds.end())
            //tell the neuron to add the synapse to this branch
            hpx_call_sync(neuronsAddr[i], DataLoader::addSynapseTarget, NULL, 0, &target, sizeof(target)) ;
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t DataLoader::addSynapseTarget = 0;
int DataLoader::addSynapseTarget_handler(const hpx_t * synapseTarget, const size_t)
{
    neurox_hpx_pin(Branch);
    assert(local->soma);
    local->soma->addSynapseTarget(*synapseTarget);
    neurox_hpx_unpin;
}

void DataLoader::registerHpxActions()
{
    neurox_hpx_register_action(1, DataLoader::createNeuron);
    neurox_hpx_register_action(1, DataLoader::addSynapseTarget);
    neurox_hpx_register_action(2, DataLoader::broadcastNetCons);
}
