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
#include "coreneuron/utils/randoms/nrnran123.h" //RNG data structures

#include "neurox/neurox.h"

using namespace std;
using namespace neurox::Input;
using namespace neurox::Input::Coreneuron;

FILE *fileNetcons;

hpx_t neuronsMutex = HPX_NULL;
std::vector<int> * neuronsGids = nullptr;

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



int DataLoader::createNeuron(NrnThread * nt, hpx_t neurond_addr)
{
        neuron_id_t neuronId = getNeuronIdFromNrnThreadId(nt->id);

        //======= 1 - reconstructs matrix (solver) values =======
        deque<Compartment*> compartments;
        for (int n=0; n<nt->end; n++)
            compartments.push_back( new Compartment(
                (offset_t) n, (floble_t) nt->_actual_a[n],
                (floble_t) nt->_actual_b[n], (floble_t) nt->_actual_d[n],
                (floble_t) nt->_actual_v[n], (floble_t) nt->_actual_rhs[n],
                (floble_t) nt->_actual_area[n],(offset_t) nt->_v_parent_index[n]));

        //reconstructs parents tree
        for (int n=1; n<nt->end; n++) //exclude top (no parent)
        {
            Compartment * parentCompartment = compartments.at(nt->_v_parent_index[n]);
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
        unsigned dataTotalOffset=nt->end*6;
        unsigned pointProcTotalOffset=0;
        map<offset_t, pair<int,offset_t>> offsetToInstance; //map of data offset -> mech instance (mech Id, node Id)
        for (NrnThreadMembList* tml = nt->tml; tml!=NULL; tml = tml->next) //For every mechanism
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
                void ** vdata = &nt->_vdata[vdataTotalOffset];
                Compartment * compartment = compartments.at(ml->nodeindices[n]);
                assert(compartment->id == ml->nodeindices[n]);
                for (int i=0; i<mech->dataSize; i++)
                {   assert(data[i]==nt->_data[dataTotalOffset+i]); }


                if (mech->pntMap > 0) //vdata
                {
                    assert((type==IClamp && mech->pdataSize==2 && mech->vdataSize==1)
                        ||((type==ProbAMPANMDA_EMS || type==ProbGABAAB_EMS) && mech->pdataSize==3 && mech->vdataSize==2));

                    //pdata[0]: offset in data (area)
                    //pdata[1]: offset for point proc in vdata[0]
                    //pdata[2]: offset for RNG in vdata[1]

                    //position vdata[0]: Point_proc
                    Point_process * ppn = &nt->pntprocs[pointProcTotalOffset++];
                    assert(nt->_vdata[pdata[1]] == ppn);
                    Point_process * pp = (Point_process*) vdata[0];
                    compartment->addSerializedVdata( (unsigned char*) (void*) pp, sizeof(Point_process));


                    //position vdata[1]: RNG (if any)
                    if (mech->vdataSize==2)
                    {
                        assert(pdata[1]+1 ==pdata[2]); //TODO no need to store offsets 1 and 2 if they are sequential
                        assert(nt->_vdata[pdata[2]] == vdata[1]);
                        nrnran123_State * RNG = (nrnran123_State*) vdata[1];
                        compartment->addSerializedVdata( (unsigned char*) (void*) RNG, sizeof(nrnran123_State));
                    }
                }
                else
                {
                    assert(mech->vdataSize==0);
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
        map< neuron_id_t, vector<NetConX*> > netcons ; //netcons per pre-synaptic neuron id)
        for (int n = 0; n < nt->n_netcon; ++n) {
            NetCon* nc = nt->netcons + n;
            assert(netcon_srcgid.size()>0); //if size==0, then setup_cleanup() in nrn_setup.cpp was called
            assert(nt->id < netcon_srcgid.size());
            assert(netcon_srcgid.at(nt->id)!=NULL);
            int srcgid = netcon_srcgid[nt->id][n];
            assert(srcgid>=0); //gids can be negative! this is a reminder that i should double-check when they happen
            int mechtype = nc->target_->_type;
            int weightscount = pnt_receive_size[mechtype];
            size_t weightsOffset = nc->u.weight_index_;
            assert(weightsOffset < nt->n_weight);


            int nodeid = nrn_threads[nc->target_->_tid]._ml_list[mechtype]->nodeindices[nc->target_->_i_instance];
            Compartment * comp = compartments.at(nodeid);
            NetConX * nx = new NetConX(mechtype, (offset_t) nc->target_->_i_instance,
                                       (floble_t) nc->delay_, weightsOffset, weightscount, nc->active_);
            double * weights = &nt->weights[weightsOffset];
            comp->addNetCon(srcgid, nx, weights);
            netcons[srcgid].push_back(nx);
        }


        if (inputParams->outputNetconsDot)
        {
          int netConsFromOthers=0;
          for (auto nc : netcons)
          {
            int srcGid = nc.first;
            if (std::find(neuronsGids->begin(), neuronsGids->end(), srcGid) == neuronsGids->end())
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


        for (auto nc : netcons)
        {
            for (auto ncv : nc.second)
                delete ncv;
            nc.second.clear();
        }
        netcons.clear();


        //======= 4 - reconstruct VecPlayContinuous events =======
        for (int v=0; v<nt->n_vecplay; v++)
        {
            VecPlayContinuous *vec = (VecPlayContinuous*) nt->_vecplay[v];
            //discover node, mechanism and data offset id that *pd points to
            PointProcInfo ppi = getPointProcInfoFromDataPointer(nt, vec->pd_);
            ppi.size = vec->y_->size();
            compartments.at(ppi.nodeId)->addVecPlay(vec->t_->data(), vec->y_->data(), ppi);
        }


        //======= 5 - recursively create branches tree ===========
        floble_t APthreshold = (floble_t) nrn_threads[nt->id].presyns[0].threshold_;
        int thvar_index = nrn_threads[nt->id].presyns[0].thvar_index_;


        createBranch(nt->id, neurond_addr, compartments, compartments.at(0),
                     (size_t) compartments.size(), offsetToInstance);

        hpx_call_sync(neurond_addr, Branch::initSoma, NULL, 0,
                      &neuronId, sizeof(neuron_id_t),
                      &APthreshold, sizeof(floble_t),
                      &thvar_index, sizeof(thvar_index));

        for (auto c : compartments)
            delete c;
        for (auto nc: netcons)
            for (auto nvc : nc.second)
                delete nvc;
}

void DataLoader::cleanCoreneuronData()
{
    nrn_cleanup();
}

void Input::Coreneuron::DataLoader::initAndLoadCoreneuronData(int argc, char ** argv)
{
    cn_input_params input_params;
    nrn_init_and_load_data(argc, argv, input_params, false);
}

int DataLoader::getMyNrnNeuronsCount()
{
    assert(nrn_threads);
    return std::accumulate(nrn_threads, nrn_threads+nrn_nthread, 0, [](int n, NrnThread & nt){return n+nt.ncell;});
}

hpx_action_t DataLoader::initMechanisms = 0;
int DataLoader::initMechanisms_handler()
{
    neurox_hpx_pin(uint64_t);

    //if serial loading and I'mnot the one who loaded the data... do nothing.
    if (!inputParams->parallelDataLoading && hpx_get_my_rank()> 0) neurox_hpx_unpin;
    int myNeuronsCount = getMyNrnNeuronsCount();
    assert(myNeuronsCount>0);

    for (int i=0; i<myNeuronsCount; i++)
    {   assert(nrn_threads[i].ncell == 1); }

    // Reconstructs unique data related to each mechanism*
    std::vector<Mechanism> mechsData;
    std::vector<int> mechsSuccessorsId, mechsDependenciesId;
    std::vector<char> mechsSym;

    //Different nrn_threads[i] have diff mechanisms sets; we'll get the union of all neurons' mechs
    std::list<NrnThreadMembList*> uniqueMechs; //list of unique mechanisms
    std::set<int> uniqueMechIds; //list of unique mechanism ids

    //insert all mechs from first neuron
    for (NrnThreadMembList* tml = nrn_threads[0].tml; tml!=NULL; tml = tml->next)
    {
        uniqueMechs.push_back(tml);
        uniqueMechIds.insert(tml->index);
    }

    //insert all mechs from other neurons that do not yet exist
    for (int i=1; i<myNeuronsCount; i++)
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
        vector<int> successorsIds, dependenciesIds;
        assert(nrn_watch_check[type] == NULL); //not supported yet

        if (inputParams->multiMex)
        {
          for (auto & tml2 : uniqueMechs)
          {
            int otherType = tml2->index;
            for (int d=0; d<nrn_prop_dparam_size_[otherType]; d++)
            {
              int ptype = memb_func[otherType].dparam_semantics[d];
              if (otherType == type && ptype>0 && ptype<1000)
                if (std::find(dependenciesIds.begin(), dependenciesIds.end(), ptype) == dependenciesIds.end())
                  dependenciesIds.push_back(ptype); //parent on dependency graph
              if (otherType!= type && ptype==type)
                if (std::find(successorsIds.begin(), successorsIds.end(), otherType) == successorsIds.end())
                  successorsIds.push_back(otherType); //children on dependency graph
            }
          }
        }
        else
        {
          //all except last have one successor
          if (tml->index != uniqueMechs.back()->index)
          {
              auto tml_next_it = std::next(tml_it,1);
              successorsIds.push_back((*tml_next_it)->index);
          }
          //all except first have one dependency
          if (tml->index != CAP)
          {
              auto tml_prev_it = std::prev(tml_it,1);
              dependenciesIds.push_back((*tml_prev_it)->index);
          }
        }

        int symLength = memb_func[type].sym ? std::strlen(memb_func[type].sym) : 0;
        if (strcmp(memb_func[type].sym, "PatternStim")==0 && inputParams->patternStim[0]=='\0')
                continue; //only load PatternStim if path initialized

        mechsData.push_back(
            Mechanism (type, nrn_prop_param_size_[type], nrn_prop_dparam_size_[type],
                       nrn_is_artificial_[type], pnt_map[type], nrn_is_ion(type),
                       symLength, nullptr, //sym will be serialized below
                       dependenciesIds.size(), nullptr, //parents will be serialized below
                       successorsIds.size(), nullptr)); //children will be serialized below

        mechsSuccessorsId.insert(mechsSuccessorsId.end(), successorsIds.begin(), successorsIds.end());
        mechsDependenciesId.insert(mechsDependenciesId.end(), dependenciesIds.begin(), dependenciesIds.end());
        mechsSym.insert(mechsSym.end(), memb_func[type].sym, memb_func[type].sym + symLength);
    }

    if (inputParams->patternStim[0]!='\0') //"initialized"
    {
        assert(0); //not an error: should be initialized already by coreneuron (and above)
        //in the future this should be contidtional (once we get rid of coreneuron data loading)
    }

    //if only one CPU (me) loaded the data, then share it with others
    if (!inputParams->parallelDataLoading && hpx_get_num_ranks()>1)
    {
      printf("neurox::setMechanisms...\n");
      vector<unsigned char> mechHasEntryInIonMap;
      vector<double> mechsMapInfo;
      for (int type=0; type<nrn_ion_global_map_size; type++)
      {
        if (nrn_ion_global_map[type]==NULL)
            mechHasEntryInIonMap.push_back(0);
        else
        {
            mechHasEntryInIonMap.push_back(1);
            mechsMapInfo.push_back(nrn_ion_global_map[type][0]);
            mechsMapInfo.push_back(nrn_ion_global_map[type][1]);
            mechsMapInfo.push_back(nrn_ion_global_map[type][2]);
        }
      }

      hpx_bcast_rsync(neurox::setMechanisms,
                    mechsData.data(), sizeof(Mechanism)*mechsData.size(),
                    mechsDependenciesId.data(), sizeof(int)* mechsDependenciesId.size(),
                    mechsSuccessorsId.data(), sizeof(int)* mechsSuccessorsId.size(),
                    mechsSym.data(), sizeof(char)*mechsSym.size());

      assert(0); //TODO need to share static _na_type_ etc inside mechanisms!x
      printf("neurox::setMechanismsGlobarVars...\n");
      hpx_bcast_rsync(neurox::setMechanismsGlobalVars,
                    &celsius, sizeof(double),
                    &nrn_ion_global_map_size, sizeof(int),
                    mechHasEntryInIonMap.data(), sizeof(unsigned char)*mechHasEntryInIonMap.size(),
                    mechsMapInfo.data(), sizeof(double)*mechsMapInfo.size());
    }
    else  //set mechanisms local to my node
    {
        setMechanisms2(mechsData.size(), mechsData.data(), mechsDependenciesId.data(),
                           mechsSuccessorsId.data(), mechsSym.data());
    }
    mechsData.clear(); mechsSuccessorsId.clear(); mechsSym.clear();

#if !defined(NDEBUG) && defined(CORENEURON_H)
    //TODO why this stopped working after HPX update? repair!
    //neurox::Input::Coreneuron::Debugger::compareMechanismsFunctionPointers(uniqueMechs);
#endif

    if (inputParams->outputMechanismsDot)
    {
      FILE *fileMechs = fopen(string("mechanisms_"+std::to_string(hpx_get_my_rank())+".dot").c_str(), "wt");
      fprintf(fileMechs, "digraph G\n{ bgcolor=%s; %s\n", DOT_PNG_BACKGROUND_COLOR,
              !inputParams->multiMex ? "layout=circo; scale=0.23;" : "");
      fprintf(fileMechs, "graph [ratio=0.3];\n", "start");
      fprintf(fileMechs, "%s [style=filled, shape=Mdiamond, fillcolor=beige];\n", "start");
      fprintf(fileMechs, "%s [style=filled, shape=Mdiamond, fillcolor=beige];\n", "end");
      fprintf(fileMechs, "\"%s (%d)\" [style=filled, fillcolor=beige];\n",
            getMechanismFromType(CAP)->sym, CAP);
      if (!inputParams->multiMex)
      {
        fprintf(fileMechs, "end -> start [color=transparent];\n");
        fprintf(fileMechs, "start -> \"%s (%d)\";\n", getMechanismFromType(CAP)->sym, CAP);
      }
      for (int m =0; m< mechanismsCount; m++)
      {
        Mechanism * mech = neurox::mechanisms[m];

        if (mech->pntMap > 0) //if is point process make it dotted
            fprintf(fileMechs, "\"%s (%d)\" [style=dashed];\n", mech->sym, mech->type);

        if (mech->dependenciesCount==0 && mech->type!=CAP) //top mechanism
            fprintf(fileMechs, "%s -> \"%s (%d)\";\n", "start", mech->sym, mech->type);

        if (mech->successorsCount==0 && mech->type!= CAP) //bottom mechanism
            fprintf(fileMechs, "\"%s (%d)\" -> %s;\n", mech->sym, mech->type, "end");

        for (int s=0; s<mech->successorsCount; s++)
        {
            Mechanism * successor = getMechanismFromType(mech->successors[s]);
            fprintf(fileMechs, "\"%s (%d)\" -> \"%s (%d)\";\n",
                    mech->sym, mech->type, successor->sym, successor->type);
        }

        if (inputParams->multiMex)
        for (int d=0; d<mech->dependenciesCount; d++)
        {
            Mechanism * parent = getMechanismFromType(mech->dependencies[d]);
            if (strcmp("SK_E2", mech->sym)==0 && strcmp("ca_ion", parent->sym)==0) continue; //TODO: hardcoded exception
            if (parent->getIonIndex() < Branch::MechanismsGraph::IonIndex::size_writeable_ions) //ie is writeable
                fprintf(fileMechs, "\"%s (%d)\" -> \"%s (%d)\" [style=dashed, arrowtype=open];\n",
                        mech->sym, mech->type, parent->sym, parent->type);
        }
      }
      fprintf(fileMechs, "}\n");
      fclose(fileMechs);
    }
    neurox_hpx_unpin;
}

hpx_action_t DataLoader::init =0;
int DataLoader::init_handler()
{
    neurox_hpx_pin(uint64_t);
    neurox::neurons = new std::vector<hpx_t>();
    neuronsGids  = new std::vector<int>();
    neuronsMutex = hpx_lco_sema_new(1);

    if (  inputParams->parallelDataLoading //disable output of netcons for parallel loading
       && inputParams->outputNetconsDot)
    {
        inputParams->outputNetconsDot=false;
        if (hpx_get_my_rank()==0)
            printf("Warning: output of netcons.dot disabled for parallel loading\n");
    }

    if (inputParams->outputNetconsDot)
    {
      assert(HPX_LOCALITIES == 1);
      fileNetcons = fopen(string("netcons.dot").c_str(), "wt");
      fprintf(fileNetcons, "digraph G\n{ bgcolor=%s; layout=circo;\n", DOT_PNG_BACKGROUND_COLOR);
      #if NETCONS_DOT_OUTPUT_NETCONS_FROM_EXTERNAL_NEURONS==true
        fprintf(fileNetcons, "external [color=gray fontcolor=gray];\n");
      #endif
    }

    neurox_hpx_unpin;
}
hpx_action_t DataLoader::initNeurons = 0;
int DataLoader::initNeurons_handler()
{
    neurox_hpx_pin(uint64_t);

    //if serial loading, and done by someone else, continue;
    if (!inputParams->parallelDataLoading && hpx_get_my_rank()> 0) neurox_hpx_unpin;

    int myNeuronsCount =  getMyNrnNeuronsCount();

#if COMPARTMENTS_DOT_OUTPUT_CORENEURON_STRUCTURE == true
    if (inputParams->outputCompartmentsDot)
    {
      for (int i=0; i<myNeuronsCount; i++)
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

    //allocate HPX memory space for neurons
    hpx_t myNeuronsGas = HPX_NULL;
    if (inputParams->parallelDataLoading) //if shared pre-balanced loading, store neurons locally
        myNeuronsGas = hpx_gas_calloc_local(myNeuronsCount, sizeof(Branch), NEUROX_HPX_MEM_ALIGNMENT);
    else //if I'm the only one loading data... spread it
        myNeuronsGas = hpx_gas_calloc_cyclic(myNeuronsCount, sizeof(Branch), NEUROX_HPX_MEM_ALIGNMENT);
    assert(myNeuronsGas != HPX_NULL);

    std::vector<int> myNeuronsGids(myNeuronsCount);
    std::vector<hpx_t> myNeurons(myNeuronsCount);
    for (int n=0; n<myNeuronsCount; n++)
    {
        myNeurons.at(n) = hpx_addr_add(myNeuronsGas, sizeof(Branch)*n, sizeof(Branch));
        myNeuronsGids.at(n)  = getNeuronIdFromNrnThreadId(n);
    }

    hpx_bcast_rsync(DataLoader::addNeurons, &myNeuronsCount, sizeof(int),
         myNeuronsGids.data(), sizeof(int)*myNeuronsGids.size(),
         myNeurons.data(), sizeof(hpx_t)*myNeurons.size());

    assert(neuronsGids->size()>0); //at least my neuron!
    for (int n=0; n<myNeuronsCount; n++)
        createNeuron(&nrn_threads[n], myNeurons.at(n));

    if (inputParams->allReduceAtLocality)
        Neuron::SlidingTimeWindow::AllReduceLocality::localityNeurons = new std::vector<hpx_t>(myNeurons);

    myNeurons.clear();
    myNeuronsGids.clear();
    neurox_hpx_unpin;
}

hpx_action_t DataLoader::addNeurons = 0;
int DataLoader::addNeurons_handler(const int nargs, const void *args[], const size_t[])
{
    /** nargs=2 where
     * args[0] = neuronsCount
     * args[1] = neurons Gids
     * args[2] = neurons hpx addr
     */

    neurox_hpx_pin(uint64_t);
    assert(nargs==3);
    const int recvNeuronsCount = *(const int*)args[0];
    const int * neuronsGids_serial = (const int*) args[1];
    const hpx_t * neuronsAddr_serial = (const hpx_t*) args[2];

    hpx_lco_sema_p(neuronsMutex);
    for (int i=0; i<recvNeuronsCount; i++)
    {
        neurons->push_back(neuronsAddr_serial[i]);
        neuronsGids->push_back(neuronsGids_serial[i]);
    }
    neurons->shrink_to_fit();
    neuronsGids->shrink_to_fit();
    hpx_lco_sema_v_sync(neuronsMutex);
    neurox_hpx_unpin;
}

hpx_action_t DataLoader::finalize = 0;
int DataLoader::finalize_handler()
{
    neurox_hpx_pin(uint64_t);

    if (inputParams->outputNetconsDot)
    {
      fprintf(fileNetcons, "}\n");
      fclose(fileNetcons);
    }

    neuronsGids->clear();  delete neuronsGids;  neuronsGids  = nullptr;
    nrn_setup_cleanup();
    hpx_lco_delete_sync(neuronsMutex);

#if defined(NDEBUG) && defined(CORENEURON_H)
    DataLoader::cleanCoreneuronData(); //if not on debug, there's no CoreNeuron comparison, so data can be cleaned-up now
#endif
    neurox_hpx_unpin;
}

void DataLoader::getNetConsBranchData(
        deque<Compartment*> & compartments, vector<NetConX> & branchNetCons,
        vector<neuron_id_t> & branchNetConsPreId, vector<floble_t> & branchWeights)
{
    for (auto & comp : compartments)
    {
        branchNetCons.insert(branchNetCons.end(), comp->netcons.begin(), comp->netcons.end());
        branchNetConsPreId.insert(branchNetConsPreId.end(), comp->netconsPreSynIds.begin(), comp->netconsPreSynIds.end());
        branchWeights.insert(branchWeights.end(), comp->netconsWeights.begin(), comp->netconsWeights.end());
    }

    //correct weighIndex variables
    int weightOffset=0;
    for (auto & nc : branchNetCons)
    {
        nc.weightIndex = weightOffset;
        weightOffset += nc.weightsCount;
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
        vector<offset_t> & pdata, vector<unsigned char> & vdata, vector<offset_t> & p,
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
    vector<vector<unsigned char>> vdataMechs (mechanismsCount);

    int n=0;
    vector<offset_t> pdataType; //type of pdata offset per pdata entry
    map< pair<int, offset_t>, offset_t> instanceToOffset; //from pair of < mech type, OLD node id> to offset in NEW representation
    map<offset_t,offset_t> fromNeuronToBranchId; //map of branch id per compartment id

    //merge all instances of all compartments into instances of the branch
    for (auto comp : compartments)
    {
        int compDataOffset=0;
        int compPdataOffset=0;
        int compVdataOffset=0;
        fromNeuronToBranchId[comp->id] = n;
        for (int m=0; m<comp->mechsTypes.size(); m++) //for all instances
        {
            int type = comp->mechsTypes[m];
            int mechOffset = mechanismsMap[type];
            assert(mechOffset>=0 && mechOffset<mechanismsCount);
            Mechanism * mech = mechanisms[mechOffset];
            dataMechs[mechOffset].insert(dataMechs[mechOffset].end(), &comp->data[compDataOffset], &comp->data[compDataOffset+mech->dataSize] );
            pdataMechs[mechOffset].insert(pdataMechs[mechOffset].end(), &comp->pdata[compPdataOffset], &comp->pdata[compPdataOffset+mech->pdataSize] );
            nodesIndicesMechs[mechOffset].push_back(comp->id);
            instancesCount[mechOffset]++;
            compDataOffset  += mech->dataSize;
            compPdataOffset += mech->pdataSize;

            if (mech->pntMap > 0) //vdata
            {
                assert((type==IClamp && mech->pdataSize==2 && mech->vdataSize==1)
                    ||((type==ProbAMPANMDA_EMS || type==ProbGABAAB_EMS) && mech->pdataSize==3 && mech->vdataSize==2));
                size_t totalVdataSize = sizeof(Point_process) + (mech->vdataSize==2 ? sizeof(nrnran123_State) : 0 );
                vdataMechs[mechOffset].insert(vdataMechs[mechOffset].end(), &comp->vdata[compVdataOffset], &comp->vdata[compVdataOffset+totalVdataSize] );
                compVdataOffset += totalVdataSize;
            }
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
            nodesIndices.push_back(fromNeuronToBranchId[ nodesIndicesMechs[m][i] ]);
            dataOffset  += mech->dataSize;
            pdataOffset += mech->pdataSize;

            if (mech->pntMap > 0) //vdata
            {
                size_t totalVdataSize = sizeof(Point_process) + (mech->vdataSize==2 ? sizeof(nrnran123_State) : 0 );
                vdata.insert(vdata.end(), &vdataMechs[m][vdataOffset], &vdataMechs[m][vdataOffset+totalVdataSize]);
                vdataOffset += totalVdataSize;
            }
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
        case -8: //"bbcorepointer"
            assert(0); //watch condition, not supported
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

hpx_t DataLoader::createBranch(int nrnThreadId, hpx_t target, deque<Compartment*> & compartments, Compartment * topCompartment,
                               int totalN, map<offset_t, pair<int,offset_t>> & offsetToInstance)
{
    assert(topCompartment!=NULL);
    offset_t n; //number of compartments in branch
    vector<floble_t> data; //compartments info (RHS, D, A, B, V, AREA)*n
    vector<offset_t> pdata; //pointers to data
    vector<offset_t> p; //parent nodes index
    vector<offset_t> instancesCount (mechanismsCount);
    vector<offset_t> nodesIndices;
    vector<hpx_t> branches;
    vector<unsigned char> vdata; //Serialized Point Processes and Random123

    //Vector Play instances
    vector<floble_t> vecPlayT;
    vector<floble_t> vecPlayY;
    vector<PointProcInfo> vecPlayInfo;

    //branch NetCons
    vector<NetConX> branchNetCons;
    vector<neuron_id_t> branchNetConsPreId;
    vector<floble_t> branchWeights;

    if (inputParams->multiSplix)
    {
        assert(0);//N/A: broken until we manage to translate mechId+instance in neuron to branch!
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
        getVecPlayBranchData(subSection, vecPlayT, vecPlayY, vecPlayInfo);
        getNetConsBranchData(subSection, branchNetCons, branchNetConsPreId, branchWeights);

        //recursively create children branches
        for (size_t c=0; c<comp->branches.size(); c++)
            branches.push_back(createBranch(nrnThreadId, HPX_NULL, compartments, comp->branches[c], totalN, offsetToInstance));
    }
    else //Flat a la Coreneuron
    {
        n = getBranchData(compartments, data, pdata, vdata, p, instancesCount, nodesIndices, totalN, offsetToInstance);
        getVecPlayBranchData(compartments, vecPlayT, vecPlayY, vecPlayInfo);
        getNetConsBranchData(compartments, branchNetCons, branchNetConsPreId, branchWeights);
    }

    //Allocate HPX Branch (top has already been created on main neurons array)
    hpx_t branchAddr = target==HPX_NULL ? hpx_gas_calloc_local(1, sizeof(Branch), NEUROX_HPX_MEM_ALIGNMENT) : target;
    bool multiSplix = inputParams->multiSplix;

    //initiate branch
    hpx_call_sync(branchAddr, Branch::init, NULL, 0,
                  &n, sizeof(offset_t),
                  &nrnThreadId, sizeof(int),
                  data.size()>0 ? data.data() : nullptr, sizeof(floble_t)*data.size(),
                  pdata.size()>0 ? pdata.data() : nullptr, sizeof(offset_t)*pdata.size(),
                  instancesCount.data(), instancesCount.size()*sizeof(offset_t),
                  nodesIndices.data(), nodesIndices.size()*sizeof(offset_t),
                  multiSplix ? branches.data() : nullptr, multiSplix ? sizeof(hpx_t)*branches.size() : 0,
                  multiSplix ? nullptr : p.data(), multiSplix ? 0 : sizeof(offset_t)*p.size(),
                  vecPlayT.size() > 0 ? vecPlayT.data() : nullptr, sizeof(floble_t)*vecPlayT.size(),
                  vecPlayY.size() > 0 ? vecPlayY.data() : nullptr, sizeof(floble_t)*vecPlayY.size(),
                  vecPlayInfo.size() > 0 ? vecPlayInfo.data() : nullptr, sizeof(PointProcInfo)*vecPlayInfo.size(),
                  branchNetCons.size() > 0 ? branchNetCons.data() : nullptr, sizeof(NetConX)*branchNetCons.size(),
                  branchNetConsPreId.size() > 0 ? branchNetConsPreId.data() : nullptr, sizeof(neuron_id_t)*branchNetConsPreId.size(),
                  branchWeights.size() > 0 ? branchWeights.data() : nullptr, sizeof(floble_t)*branchWeights.size(),
                  vdata.size()>0 ? vdata.data() : nullptr, sizeof(unsigned char)*vdata.size()
                  );

    return branchAddr;
}

hpx_action_t DataLoader::initNetcons = 0;
int DataLoader::initNetcons_handler()
{
    neurox_hpx_pin(Branch);
//    neurox_hpx_recursive_branch_async_call(DataLoader::initNetcons);

    hpx_addr_t lco_branches = local->branchTree && local->branchTree->branchesCount ? hpx_lco_and_new(local->branchTree->branchesCount) : HPX_NULL;
    if (local->branchTree)
    for (int c=0; c<local->branchTree->branchesCount; c++)
       {_hpx_call(local->branchTree->branches[c], DataLoader::initNetcons, lco_branches, NULL, 0);}


    if (inputParams->outputNetconsDot)
        fprintf(fileNetcons, "%d [style=filled, shape=ellipse];\n", local->soma->gid);

    const floble_t impossiblyLargeDelay = 99999999;
    for (int i=0; i < neurox::neurons->size(); i++) //loop through all neurons
    {
        neuron_id_t srcGid = neuronsGids->at(i);

        //if I'm connected to it (ie is not artificial or non-existent)
        if (local->netcons.find(srcGid) != local->netcons.end())
        {
            hpx_t srcAddr = neurox::neurons->at(i);
            floble_t minDelay = impossiblyLargeDelay;
            for (NetConX *& nc : local->netcons.at(srcGid))
                if (nc->active)
                   minDelay = min(minDelay, nc->delay);

#if NETCONS_OUTPUT_ADDITIONAL_VALIDATION_FILE==true
            if (inputParams->outputNetconsDot && minDelay!=impossiblyLargeDelay)
                fprintf(fileNetcons, "%d -> %d [label=\"%d (%.2fms)\"];\n",
                    srcGid, local->soma->gid, local->netcons.at(srcGid).size(), minDelay);
#endif

            //tell the neuron to add a synapse to this branch
            //and inform him of the fastest netcon we have
            if (minDelay != impossiblyLargeDelay) //if any active synapse
            {
                //inform pre-syn neuron that he connects to me
                hpx_call_sync(srcAddr, DataLoader::addSynapse, NULL, 0,
                    &target, sizeof(target), &minDelay, sizeof(minDelay),
                    &local->soma->gid, sizeof(int));

                //add this pre-syn neuron as my time-dependency
                assert(local->soma);
                if (inputParams->algorithm == Algorithm::BackwardEulerWithPairwiseSteping)
                  local->soma->timeDependencies->updateTimeDependency(srcGid, local->soma->gid,
                       inputParams->tstart+minDelay*Neuron::TimeDependencies::notificationIntervalRatio, true);
            }
        }
    }
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t DataLoader::addSynapse = 0;
int DataLoader::addSynapse_handler(
        const int nargs, const void *args[], const size_t[] )
{
    neurox_hpx_pin(Branch);
    assert(local->soma);
    assert(nargs==3);
    hpx_t addr = *(const hpx_t*) args[0];
    floble_t minDelay = *(const floble_t*) args[1];
    int destinationGid = *(const int*) args[2];
    local->soma->addSynapse(new Neuron::Synapse(addr,minDelay,destinationGid));
    neurox_hpx_unpin;
}

void DataLoader::registerHpxActions()
{
    neurox_hpx_register_action(0, DataLoader::init);
    neurox_hpx_register_action(0, DataLoader::initMechanisms);
    neurox_hpx_register_action(0, DataLoader::initNeurons);
    neurox_hpx_register_action(0, DataLoader::initNetcons);
    neurox_hpx_register_action(0, DataLoader::finalize);
    neurox_hpx_register_action(2, DataLoader::addSynapse);
    neurox_hpx_register_action(2, DataLoader::addNeurons);
    neurox_hpx_register_action(0, DataLoader::initNetcons);
}