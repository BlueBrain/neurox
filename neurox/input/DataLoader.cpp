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
using namespace neurox::input;
using namespace neurox::algorithms;

static FILE *fileNetcons;
static hpx_t neuronsMutex = HPX_NULL;
static std::vector<int> * neuronsGids = nullptr;

static tools::LoadBalancing* loadBalancing = nullptr;

neuron_id_t DataLoader::GetNeuronIdFromNrnThreadId(int nrn_id)
{
    return (neuron_id_t) nrn_threads[nrn_id].presyns[0].gid_;
}

void DataLoader::PrintSubClustersToFile(FILE * fileCompartments, Compartment *topCompartment)
{
  if (inputParams->outputCompartmentsDot)
  {
    assert(topCompartment!=NULL);
    fprintf(fileCompartments, "subgraph cluster_%d { ", topCompartment->id);
    if (topCompartment->id == 0)
        fprintf(fileCompartments, "label=\"SOMA\"; ");
    else if (topCompartment->id == 2 )
        fprintf(fileCompartments, "label=\"AIS\"; ");
    Compartment * comp=NULL;
    for (comp = topCompartment; comp->branches.size()==1; comp = comp->branches.front())
        fprintf(fileCompartments, "%d; ", comp->id);
    fprintf(fileCompartments, "%d };\n", comp->id);
    for (int c=0; c<comp->branches.size(); c++)
        PrintSubClustersToFile(fileCompartments, comp->branches[c]);
  }
}

PointProcInfo DataLoader::GetPointProcInfoFromDataPointer(NrnThread * nt, double *pd, size_t size)
{
    PointProcInfo ppi;
    ppi.nodeId=-1;
    ppi.size = size;
    bool found=false;
    for (NrnThreadMembList* tml = nt->tml; !found && tml!=NULL; tml = tml->next) //For every mechanism
    {
        int type = tml->index;
        Mechanism * mech = GetMechanismFromType(type);
        Memb_list * ml = tml->ml;
        for (int n=0; n<ml->nodecount; n++)
            for (int i=0; i<mech->dataSize; i++)
            {
#if LAYOUT==1
                int dataOffset = mech->dataSize*n+i;
#else
                int dataOffset = tools::Vectorizer::SizeOf(ml->nodecount)*i+n ;
#endif
                if (&ml->data[dataOffset] != pd) continue; //if not this variable

                ppi.nodeId = ml->nodeindices[n];
                ppi.mechType = type;
                ppi.mechInstance = (offset_t) n;
                ppi.instanceDataOffset = i;
                found=true;
                break;
           }
    }
    assert(found);
    return ppi;
}



int DataLoader::CreateNeuron(int neuron_idx, void * targets)
{
        NrnThread * nt = &nrn_threads[neuron_idx];
        hpx_t neuron_addr = ((hpx_t*)targets)[neuron_idx];
        neuron_id_t neuronId = GetNeuronIdFromNrnThreadId(nt->id);
        int N = nt->end;

        //if data is permuted, this method fails.
        assert(!use_interleave_permute && !use_solve_interleave && nt->_permute==NULL);

        //map of padded to non-padded offsets of data
        size_t dataSizePadded  = 6*tools::Vectorizer::SizeOf(N);
        for (NrnThreadMembList* tml = nt->tml; tml!=NULL; tml = tml->next)
            dataSizePadded += tools::Vectorizer::SizeOf(tml->ml->nodecount) * mechanisms[mechanismsMap[tml->index]]->dataSize;

        //map of post- to pre-padding values of pdata
        //only used for branched neurons, otherwise pointers for padded and non-padded layouts are the same
        std::vector<int> dataOffsets(inputParams->branchingDepth>0 ? dataSizePadded : 0, -99999);

        if (inputParams->branchingDepth>0)
          for (int n=0; n<N; n++)
            for (int i=0; i<6; i++)
            {
                int offsetPadded = tools::Vectorizer::SizeOf(N)*i+n;
                int offsetNonPadded = N*i+n;
                dataOffsets[offsetPadded] = offsetNonPadded;
              }

        //======= 1 - reconstructs matrix (solver) values =======
        deque<Compartment*> compartments;
        for (int n=0; n<nt->end; n++)
            compartments.push_back( new Compartment(
                (offset_t) n, (floble_t) nt->_actual_a[n],
                (floble_t) nt->_actual_b[n], (floble_t) nt->_actual_d[n],
                (floble_t) nt->_actual_v[n], (floble_t) nt->_actual_rhs[n],
                (floble_t) nt->_actual_area[n],(offset_t) nt->_v_parent_index[n]));

        //reconstructs parents tree
        for (int n=1; n<N; n++) //exclude top (no parent)
        {
            Compartment * parentCompartment = compartments.at(nt->_v_parent_index[n]);
            parentCompartment->AddChild(compartments.at(n));
        }

        if (inputParams->outputCompartmentsDot)
        {
          set<int> axonInitSegmentCompartments;
          FILE *fileCompartments = fopen(string("compartments_"+to_string(neuronId)+".dot").c_str(), "wt");
          fprintf(fileCompartments, "graph G_%d\n{ bgcolor=%s; \n", neuronId, DOT_PNG_BACKGROUND_COLOR );
          fprintf(fileCompartments, "graph [fontname=helvetica, style=filled, color=blue, fillcolor=floralwhite];\n");
          fprintf(fileCompartments, "node [fontname=helvetica, shape=cylinder, color=gray, style=filled, fillcolor=white];\n");
          fprintf(fileCompartments, "edge [fontname=helvetica, color=gray];\n");
          PrintSubClustersToFile(fileCompartments, compartments.at(0)); //add subclusters
          for (Compartment *comp : compartments) //draw edges
            for (int c=0; c<comp->branches.size(); c++)
            {
                bool isSoma = comp->id == 1; //bottom of soma
                Compartment * child = comp->branches.at(c);
                if ((isSoma && c==0) //connection from some to AIS
                ||  axonInitSegmentCompartments.find(comp->id) != axonInitSegmentCompartments.end()) //connection to any AIS compartment
                {
                    //reverse connection, so that it plots AIS on top of soma in dot file
                    fprintf(fileCompartments, "%d -- %d%s;\n", child->id, comp->id, comp->branches.size()==1 ? "" : " [color=black]");
                    axonInitSegmentCompartments.insert(child->id);
                }
                else
                    fprintf(fileCompartments, "%d -- %d%s;\n", comp->id, child->id, comp->branches.size()==1 ? "" : " [color=black]");
            }
          fprintf(fileCompartments, "}\n");
          fclose(fileCompartments);
        }

        //======= 2 - reconstructs mechanisms instances ========
        unsigned vdataTotalOffset=0;
        unsigned dataTotalOffset = N*6; //no padding
        unsigned dataTotalPaddedOffset = tools::Vectorizer::SizeOf(N)*6; //with padding
        unsigned pointProcTotalOffset=0;

        //information about offsets in data and node ifs of all instances of all ions
        vector<DataLoader::IonInstancesInfo> ionsInstancesInfo (Mechanism::Ion::size_all_ions);
        for (NrnThreadMembList* tml = nt->tml; tml!=NULL; tml = tml->next) //For every mechanism
        {
            int type = tml->index;
            Memb_list * ml = tml->ml; //Mechanisms application to each compartment
            Mechanism * mech = GetMechanismFromType(type);
            int ionOffset = mech->GetIonIndex();
            for (int n=0; n<ml->nodecount; n++) //for every mech instance (or compartment this mech is applied to)
            {
                if (mech->isIon)
                {
                    if (n==0)
                    {
                        ionsInstancesInfo[ionOffset].mechType  = type;
                        ionsInstancesInfo[ionOffset].dataStart = dataTotalOffset;
                        ionsInstancesInfo[ionOffset].dataEnd   = dataTotalOffset + ml->nodecount*mech->dataSize;
                    }
                    ionsInstancesInfo[ionOffset].nodeIds.push_back(ml->nodeindices[n]);
                }

                std::vector<double> data;
                for (int i=0; i<mech->dataSize; i++)
                {
                    int offsetNonPadded = mech->dataSize*n+i;
#if LAYOUT==1
                    data.push_back(ml->data[offsetNonPadded]);
                    assert(ml->data[offsetNonPadded]==nt->_data[dataTotalOffset+offsetNonPadded]);
#else
                    int offsetPadded = tools::Vectorizer::SizeOf(ml->nodecount)*i+n ;
                    data.push_back(ml->data[offsetPadded]);
                    assert(ml->data[offsetPadded]==nt->_data[dataTotalPaddedOffset+offsetPadded]);
                    if (inputParams->branchingDepth>0)
                        dataOffsets[dataTotalPaddedOffset+offsetPadded] = dataTotalOffset+offsetNonPadded;
#endif
                }
                data.shrink_to_fit();

                std::vector<int> pdata;
                for (int i=0; i<mech->pdataSize; i++)
                {
#if LAYOUT==1
                    int pdataOffsetNonPadded = mech->pdataSize*n+i;
                    pdata.push_back(ml->pdata[pdataOffsetNonPadded]);
#else
                    int pdataOffsetPadded = tools::Vectorizer::SizeOf(ml->nodecount)*i+n ;
                    int pd = ml->pdata[pdataOffsetPadded];
                    int ptype = memb_func[mech->type].dparam_semantics[i];

                    //remove extra space added by padding (for pointer to area or ion mech instance)
                    if (inputParams->branchingDepth>0 && (ptype==-1 || (ptype>0 && ptype<1000)))
                    {
                        assert(dataOffsets.at(pd)!=-99999);
                        pdata.push_back(dataOffsets.at(pd)); //offset to non-padded SoA value
                    }
                    else
                        pdata.push_back(pd);
#endif
                }
                pdata.shrink_to_fit();

                void ** vdata = &nt->_vdata[vdataTotalOffset];
                Compartment * compartment = compartments.at(ml->nodeindices[n]);
                assert(compartment->id == ml->nodeindices[n]);

                if (mech->pntMap > 0 || mech->vdataSize>0) //vdata
                {
                    assert(( type == IClamp && mech->vdataSize == 1 && mech->pdataSize == 2 && mech->pntMap>0)
                        || ( type == StochKv && mech->vdataSize == 1 && mech->pdataSize == 5 && mech->pntMap==0)
                        || ((type == ProbAMPANMDA_EMS || type == ProbGABAAB_EMS)
                            && mech->vdataSize == 2 && mech->pdataSize == 3 && mech->pntMap>0 ));

                    //ProbAMPANMDA_EMS, ProbAMPANMDA_EMS and IClamp:
                    //pdata[0]: offset in data (area)
                    //pdata[1]: offset for Point_process in vdata[0]
                    //pdata[2]: offset for RNG in vdata[1]   (NOT for IClamp,  =pdata[1]+1)

                    //StochKv:
                    //pdata[0]: offset in area (ion_ek)
                    //pdata[1]: offset in area (ion_ik)
                    //pdata[2]: offset in area (ion_dikdv)
                    //pdata[3]: offset for RNG in vdata[0]
                    //pdata[4]: offset in data (area)

                    if (type == IClamp || type == ProbAMPANMDA_EMS || type == ProbGABAAB_EMS)
                    {
                        const int pointProcOffsetInPdata = 1;
                        Point_process * ppn = &nt->pntprocs[pointProcTotalOffset++];
                        assert(nt->_vdata[pdata[pointProcOffsetInPdata]] == ppn);
                        Point_process * pp = (Point_process*) vdata[0];
                        assert(pp!=NULL);
                        compartment->AddSerializedVdata( (unsigned char*) (void*) pp, sizeof(Point_process));
                    }

                    if (type == StochKv || type == ProbAMPANMDA_EMS || type == ProbGABAAB_EMS)
                    {
                        int rngOffsetInPdata = mech->type == StochKv ? 3 : 2;
                        int rngOffsetInVdata = mech->type == StochKv ? 0 : 1;
                        assert(nt->_vdata[pdata[rngOffsetInPdata]] == vdata[rngOffsetInVdata]);
                        nrnran123_State * RNG = (nrnran123_State*) vdata[rngOffsetInVdata];

                        //TODO: manual hack: StochKv's current state has NULL pointers, why?
                        if (RNG==NULL && type == StochKv)
                            compartment->AddSerializedVdata( (unsigned char*) new nrnran123_State, sizeof(nrnran123_State));
                        else
                        {
                            assert(RNG!=NULL);
                            compartment->AddSerializedVdata( (unsigned char*) (void*) RNG, sizeof(nrnran123_State));
                        }
                    }
                }
                compartment->AddMechanismInstance(type, n, data.data(), mech->dataSize, pdata.data(),  mech->pdataSize);

                vdataTotalOffset += (unsigned) mech->vdataSize;
            }
            dataTotalOffset       += mech->dataSize*ml->nodecount;
            dataTotalPaddedOffset += mech->dataSize*tools::Vectorizer::SizeOf(ml->nodecount);
        }
        for (Compartment* comp : compartments)  comp->ShrinkToFit();

        dataOffsets.clear();

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
            comp->AddNetCon(srcgid, nx, weights);
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
        for (Compartment* comp : compartments)  comp->ShrinkToFit();

        //======= 4 - reconstruct VecPlayContinuous events =======
        for (int v=0; v<nt->n_vecplay; v++)
        {
            VecPlayContinuous *vec = (VecPlayContinuous*) nt->_vecplay[v];
            //discover node, mechanism and data offset id that *pd points to
            PointProcInfo ppi = GetPointProcInfoFromDataPointer(nt, vec->pd_,  vec->y_->size());
            compartments.at(ppi.nodeId)->AddVecPlay(vec->t_->data(), vec->y_->data(), ppi);
        }


        //======= 5 - recursively create branches tree ===========
        floble_t APthreshold = (floble_t) nrn_threads[nt->id].presyns[0].threshold_;
        int thvar_index = nrn_threads[nt->id].presyns[0].thvar_index_;

        for (Compartment* comp : compartments)
            comp->ShrinkToFit();

        CreateBranch(nt->id, neuron_addr, BranchType::Soma,
                     thvar_index, compartments, compartments.at(0),
                     ionsInstancesInfo, inputParams->branchingDepth);

        hpx_call_sync(neuron_addr, Branch::InitSoma, NULL, 0,
                      &neuronId, sizeof(neuron_id_t),
                      &APthreshold, sizeof(floble_t));

        for (auto c : compartments)
            delete c;
        for (auto nc: netcons)
            for (auto nvc : nc.second)
                delete nvc;
    return 0;
}

void DataLoader::CleanCoreneuronData(const bool clean_ion_global_map)
{
    nrn_cleanup(clean_ion_global_map);
}

void DataLoader::InitAndLoadCoreneuronData(
        int argc, char ** argv,
        bool nrnmpi_under_nrncontrol,
        bool run_setup_cleanup)
{
    //nrnmpi_under_nrncontrol=true allows parallel data loading without "-m" flag, see rnmpi_init()
    nrn_init_and_load_data(argc, argv, nrnmpi_under_nrncontrol, run_setup_cleanup);
}

int DataLoader::GetMyNrnNeuronsCount()
{
    assert(nrn_threads);
    return std::accumulate(nrn_threads, nrn_threads+nrn_nthread, 0, [](int n, NrnThread & nt){return n+nt.ncell;});
}

hpx_action_t DataLoader::InitMechanisms = 0;
int DataLoader::InitMechanisms_handler()
{
    neurox_hpx_pin(uint64_t);

    //Build Mechanisms data (without dependencies) from memb_func
    neurox::mechanismsCount=0;
    neurox::mechanismsMap = new int[n_memb_func];
    for (int type=0; type<n_memb_func; type++)
    {
        neurox::mechanismsMap[type]=-1;
        if (memb_func[type].alloc == NULL) continue;
        neurox::mechanismsMap[type] = neurox::mechanismsCount++;
    }

    neurox::mechanisms = new Mechanism*[neurox::mechanismsCount];
    for (int type=0; type<n_memb_func; type++)
    {
        if (neurox::mechanismsMap[type]==-1) continue;

        int symLength = memb_func[type].sym ? std::strlen(memb_func[type].sym) : 0;
        mechanisms[neurox::mechanismsMap[type]] = new
            Mechanism (type, nrn_prop_param_size_[type], nrn_prop_dparam_size_[type],
                       nrn_is_artificial_[type], pnt_map[type], nrn_is_ion(type),
                       symLength, memb_func[type].sym, memb_func[type],
                       0, nullptr,  //parents will be communicated below
                       0, nullptr); //children will be communicated below
    }

    //build dependencies between mechanisms (only available from the morphology)
    int myNeuronsCount = GetMyNrnNeuronsCount();

    if (myNeuronsCount==0) neurox_hpx_unpin;

    for (int i=0; i<myNeuronsCount; i++)
    {   assert(nrn_threads[i].ncell == 1); }

    // Reconstructs unique data related to each mechanism
    std::vector<int> successorsCount(neurox::mechanismsCount,0), dependenciesCount(neurox::mechanismsCount,0);
    std::vector<int> successorsIds, dependenciesIds;

    //Different nrn_threads[i] have diff mechanisms; we'll get the union of all neurons' mechs
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
        vector<int> successors, dependencies;
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
                if (std::find(dependencies.begin(), dependencies.end(), ptype) == dependencies.end())
                  dependencies.push_back(ptype); //parent on dependency graph
              if (otherType!= type && ptype==type)
                if (std::find(successors.begin(), successors.end(), otherType) == successors.end())
                  successors.push_back(otherType); //children on dependency graph
            }
          }
        }
        else
        {
          //all except last have one successor
          if (tml->index != uniqueMechs.back()->index)
          {
              auto tml_next_it = std::next(tml_it,1);
              successors.push_back((*tml_next_it)->index);
          }
          //all except first have one dependency
          if (tml->index != CAP)
          {
              auto tml_prev_it = std::prev(tml_it,1);
              dependencies.push_back((*tml_prev_it)->index);
          }
        }

        dependenciesCount.at(neurox::mechanismsMap[type]) = dependencies.size();
        dependenciesIds.insert(dependenciesIds.end(), dependencies.begin(), dependencies.end());

        successorsCount.at(neurox::mechanismsMap[type]) = successors.size();
        successorsIds.insert(successorsIds.end(), successors.begin(), successors.end());
    }

    if (inputParams->patternStim[0]!='\0') //"initialized"
    {
        assert(0); //not an error: should be initialized already by coreneuron (and above)
        //in the future this should be contidtional (once we get rid of coreneuron data loading)
    }

    //set mechanisms dependencies
    if (neurox::ParallelExecution() && inputParams->branchingDepth>0)
    {
        //broadcast dependencies, most complete dependency graph will be used across the network
        //(this solves issue of localities loading morphologies without all mechanisms,
        //and processing branches of other localities with the missing mechanisms)
        hpx_bcast_rsync(DataLoader::UpdateMechanismsDependencies,
                        dependenciesCount.data(), sizeof(int)*dependenciesCount.size(),
                        dependenciesIds.data(),   sizeof(int)*dependenciesIds.size(),
                        successorsCount.data(),   sizeof(int)*successorsCount.size(),
                        successorsIds.data(),     sizeof(int)*successorsIds.size()
                        );
    }
    else
    {
        //regular setting of mechanisms: all localities have the dependency graph for their morphologies
        neurox::SetMechanismsDependencies(
                dependenciesCount.data(), dependenciesIds.data(),
                successorsCount.data(),   successorsIds.data());
    }
    neurox_hpx_unpin;
}

hpx_action_t DataLoader::Init =0;
int DataLoader::Init_handler()
{
    neurox_hpx_pin(uint64_t);
    neurox::neurons = new std::vector<hpx_t>();
    neuronsGids  = new std::vector<int>();
    neuronsMutex = hpx_lco_sema_new(1);

    if (hpx_get_my_rank()==0 && inputParams->branchingDepth>0)
        loadBalancing = new tools::LoadBalancing();

    if (neurox::ParallelExecution() //disable output of netcons for parallel loading
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

    if (inputParams->outputMechanismsDot)
    {
      FILE *fileMechs = fopen(string("mechanisms_"+std::to_string(hpx_get_my_rank())+".dot").c_str(), "wt");
      fprintf(fileMechs, "digraph G\n{ bgcolor=%s; %s\n", DOT_PNG_BACKGROUND_COLOR,
              !inputParams->multiMex ? "layout=circo; scale=0.23;" : "");
      fprintf(fileMechs, "graph [ratio=0.3];\n", "start");
      fprintf(fileMechs, "%s [style=filled, shape=Mdiamond, fillcolor=beige];\n", "start");
      fprintf(fileMechs, "%s [style=filled, shape=Mdiamond, fillcolor=beige];\n", "end");
      fprintf(fileMechs, "\"%s (%d)\" [style=filled, fillcolor=beige];\n",
            GetMechanismFromType(CAP)->membFunc.sym, CAP);
      if (!inputParams->multiMex)
      {
        fprintf(fileMechs, "end -> start [color=transparent];\n");
        fprintf(fileMechs, "start -> \"%s (%d)\";\n", GetMechanismFromType(CAP)->membFunc.sym, CAP);
      }
      for (int m =0; m< mechanismsCount; m++)
      {
        Mechanism * mech = neurox::mechanisms[m];

        if (mech->pntMap > 0) //if is point process make it dotted
            fprintf(fileMechs, "\"%s (%d)\" [style=dashed];\n", mech->membFunc.sym, mech->type);

        if (mech->dependenciesCount==0 && mech->type!=CAP) //top mechanism
            fprintf(fileMechs, "%s -> \"%s (%d)\";\n", "start", mech->membFunc.sym, mech->type);

        if (mech->successorsCount==0 && mech->type!= CAP) //bottom mechanism
            fprintf(fileMechs, "\"%s (%d)\" -> %s;\n", mech->membFunc.sym, mech->type, "end");

        for (int s=0; s<mech->successorsCount; s++)
        {
            Mechanism * successor = GetMechanismFromType(mech->successors[s]);
            fprintf(fileMechs, "\"%s (%d)\" -> \"%s (%d)\";\n",
                    mech->membFunc.sym, mech->type, successor->membFunc.sym, successor->type);
        }

        if (inputParams->multiMex)
        for (int d=0; d<mech->dependenciesCount; d++)
        {
            Mechanism * parent = GetMechanismFromType(mech->dependencies[d]);
            if (strcmp("SK_E2", mech->membFunc.sym)==0 && strcmp("ca_ion", parent->membFunc.sym)==0) continue; //TODO: hardcoded exception
            if (parent->GetIonIndex() < Mechanism::Ion::size_writeable_ions) //ie is writeable
                fprintf(fileMechs, "\"%s (%d)\" -> \"%s (%d)\" [style=dashed, arrowtype=open];\n",
                        mech->membFunc.sym, mech->type, parent->membFunc.sym, parent->type);
        }
      }
      fprintf(fileMechs, "}\n");
      fclose(fileMechs);
    }

    neurox_hpx_unpin;
}

hpx_action_t DataLoader::InitNeurons = 0;
int DataLoader::InitNeurons_handler()
{
    neurox_hpx_pin(uint64_t);

    int myNeuronsCount =  GetMyNrnNeuronsCount();

    if (myNeuronsCount==0) neurox_hpx_unpin;

#if COMPARTMENTS_DOT_OUTPUT_CORENEURON_STRUCTURE == true
    if (inputParams->outputCompartmentsDot)
    {
      for (int i=0; i<myNeuronsCount; i++)
      {
        neuron_id_t neuronId = GetNeuronIdFromNrnThreadId(i);
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
    if (neurox::ParallelExecution()) //if shared pre-balanced loading, store neurons locally
        myNeuronsGas = hpx_gas_calloc_local(myNeuronsCount, sizeof(Branch), NEUROX_MEM_ALIGNMENT);
    else //if I'm the only one loading data... spread it
        myNeuronsGas = hpx_gas_calloc_cyclic(myNeuronsCount, sizeof(Branch), NEUROX_MEM_ALIGNMENT);
    assert(myNeuronsGas != HPX_NULL);

    int * myNeuronsGids = new int[myNeuronsCount];
    hpx_t * myNeuronsTargets = new hpx_t[myNeuronsCount];

    for (int n=0; n<myNeuronsCount; n++)
    {
        myNeuronsTargets[n] = hpx_addr_add(myNeuronsGas, sizeof(Branch)*n, sizeof(Branch));
        myNeuronsGids[n] = GetNeuronIdFromNrnThreadId(n);
    }

    hpx_bcast_rsync(DataLoader::AddNeurons, &myNeuronsCount, sizeof(int),
         myNeuronsGids, sizeof(int)*myNeuronsCount,
         myNeuronsTargets, sizeof(hpx_t)*myNeuronsCount);

    //If there is branch parallelism, neurons are generated in serial
    //(otherwise branches benchmark is affected by scheduler and running threads)
    //TODO I dont think so, a benchmark thread will run until it finishes!
    if (inputParams->branchingDepth>0)
        for (int i=0; i<myNeuronsCount; i++)
            CreateNeuron(i, myNeuronsTargets);
    else
        hpx_par_for_sync( DataLoader::CreateNeuron, 0, myNeuronsCount, myNeuronsTargets);

    if (inputParams->allReduceAtLocality)
        AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::localityNeurons =
            new std::vector<hpx_t>(myNeuronsTargets, myNeuronsTargets+myNeuronsCount);

    delete [] myNeuronsGids;
    delete [] myNeuronsTargets;
    neurox_hpx_unpin;
}

hpx_action_t DataLoader::AddNeurons = 0;
int DataLoader::AddNeurons_handler(const int nargs, const void *args[], const size_t[])
{
    /**
     * nargs=3 where
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

static hpx_t dependenciesMutex = HPX_NULL;
static int maxDependenciesCount = 0;
static int maxSuccessorsCount = 0;

hpx_action_t DataLoader::UpdateMechanismsDependencies = 0;
int DataLoader::UpdateMechanismsDependencies_handler(const int nargs, const void *args[], const size_t[])
{
    /**
     * nargs=4 where
     * args[0] = successors count per mechanism
     * args[1] = successors ids
     * args[2] = dependendencies count per mechanism
     * args[3] = dependendencies ids
     */

    neurox_hpx_pin(uint64_t);
    assert(nargs==4);

    const int * dependenciesCount = (const int*)args[0];
    const int * dependenciesIds   = (const int*)args[1];
    const int * successorsCount   = (const int*)args[2];
    const int * successorsIds     = (const int*)args[3];

    hpx_lco_sema_p(neuronsMutex);
    int dependenciesTotalCount=0, successorsTotalCount=0;

    for (int m=0; m<neurox::mechanismsCount; m++)
    {
        dependenciesTotalCount += dependenciesCount[m];
        successorsTotalCount   += successorsCount[m];
    };

    //if there are more dependencies that the ones I know
    if (dependenciesTotalCount>maxDependenciesCount)
    {
        maxDependenciesCount = dependenciesTotalCount;
        neurox::SetMechanismsDependencies(dependenciesCount, dependenciesIds, nullptr, nullptr);
    }

    if (successorsTotalCount>maxSuccessorsCount)
    {
        maxSuccessorsCount = successorsTotalCount;
        neurox::SetMechanismsDependencies(nullptr, nullptr, successorsCount, successorsIds);
    }
    hpx_lco_sema_v_sync(neuronsMutex);

    neurox_hpx_unpin;
}
hpx_action_t DataLoader::Finalize = 0;
int DataLoader::Finalize_handler()
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

#if defined(NDEBUG) 
    DataLoader::CleanCoreneuronData(); //if not on debug, there's no CoreNeuron comparison, so data can be cleaned-up now
#else
    //print Load Balancing table
    if (hpx_get_my_rank()==0 && inputParams->branchingDepth>0)
        loadBalancing->PrintTable();
#endif
    delete loadBalancing; loadBalancing = nullptr;
    neurox_hpx_unpin;
}


void DataLoader::GetAllChildrenCompartments(deque<Compartment*> & subSection, Compartment * topCompartment)
{
    //parent added first, to respect solvers logic, of parents' ids first
    subSection.push_back(topCompartment);
    for (int c=0; c<topCompartment->branches.size(); c++)
        GetAllChildrenCompartments(subSection, topCompartment->branches.at(c));
}

void DataLoader::GetMechInstanceMap(deque<Compartment*> & compartments, vector<map<int,int>> & mechsInstancesMap)
{
    vector<deque<int>> mechsInstancesIds(neurox::mechanismsCount); //mech offset -> list of mech instance id
    for (Compartment * comp : compartments)
        for (int m=0; m<comp->mechsTypes.size(); m++) //for all instances
        {
            int type = comp->mechsTypes.at(m);
            int mechOffset = neurox::mechanismsMap[type];
            mechsInstancesIds[mechOffset].push_back(comp->mechsInstances[m]);
        }

    //convert neuron mech-instances ids from neuron- to branch-level
    for (int m=0; m<neurox::mechanismsCount; m++)
      for (int i=0; i<mechsInstancesIds[m].size(); i++)
      {
        int oldInstanceId = mechsInstancesIds.at(m).at(i);
        mechsInstancesMap[m][oldInstanceId] = i;
      }
}

void DataLoader::GetNetConsBranchData(
        deque<Compartment*> & compartments, vector<NetConX> & branchNetCons,
        vector<neuron_id_t> & branchNetConsPreId, vector<floble_t> & branchWeights,
        vector<map<int,int>> * mechInstanceMap)
{
    for (auto & comp : compartments)
    {
        branchNetCons.insert(branchNetCons.end(), comp->netcons.begin(), comp->netcons.end());
        branchNetConsPreId.insert(branchNetConsPreId.end(), comp->netconsPreSynIds.begin(), comp->netconsPreSynIds.end());
        branchWeights.insert(branchWeights.end(), comp->netconsWeights.begin(), comp->netconsWeights.end());
    }

    //correct weighIndex variables
    int weightOffset=0;
    for (NetConX & nc : branchNetCons)
    {
        nc.weightIndex = weightOffset;
        weightOffset += nc.weightsCount;
    }

    //convert mech instance id from neuron to branch level
    if (mechInstanceMap)
        for (NetConX & nc : branchNetCons)
            nc.mechInstance = (*mechInstanceMap)[neurox::mechanismsMap[nc.mechType]][nc.mechInstance];
}

void DataLoader::GetVecPlayBranchData(
        deque<Compartment*> & compartments, vector<floble_t> & vecPlayTdata,
        vector<floble_t> & vecPlayYdata, vector<PointProcInfo> & vecPlayInfo,
        vector<map<int,int>> * mechInstanceMap)
{
    //convert node id and mech instance id in PointProcess from neuron to branch level
    if (mechInstanceMap)
    {
        std::map<int,int> fromOldToNewCompartmentId;
        for (int n=0; n<compartments.size(); n++)
            fromOldToNewCompartmentId[compartments.at(n)->id] = n;

        for (int p=0; p<vecPlayInfo.size(); p++)
        {
            PointProcInfo & ppi = vecPlayInfo[p];
            ppi.mechInstance = (offset_t) (*mechInstanceMap)[neurox::mechanismsMap[ppi.mechType]][ppi.mechInstance];
            ppi.nodeId = fromOldToNewCompartmentId[ppi.nodeId];
        }
    }

    for (auto comp : compartments)
    {
        vecPlayTdata.insert(vecPlayTdata.end(), comp->vecPlayTdata.begin(), comp->vecPlayTdata.end());
        vecPlayYdata.insert(vecPlayYdata.end(), comp->vecPlayYdata.begin(), comp->vecPlayYdata.end());
        vecPlayInfo.insert(vecPlayInfo.end(), comp->vecPlayInfo.begin(), comp->vecPlayInfo.end());
    }
}

int DataLoader::GetBranchData(
        deque<Compartment*> & compartments, vector<floble_t> & data,
        vector<offset_t> & pdata, vector<unsigned char> & vdata, vector<offset_t> & p,
        vector<offset_t> & instancesCount, vector<offset_t> & nodesIndices,
        int N,  vector<DataLoader::IonInstancesInfo> & ionsInstancesInfo,
        vector<map<int,int>> * mechInstanceMap)
{
    for (auto comp : compartments)
        { assert (comp != NULL); }

    int n=0; //number of compartments
    int vdataPointerOffset=0;

    ////// Basic information for RHS, D, A, B, V and area
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
    for (auto comp : compartments)
        p.push_back(comp->p);

    ////// Tree of neurons: convert from neuron- to branch-level
    std::map<int,int> fromOldToNewCompartmentId;
    for (Compartment * comp : compartments)
    {
        fromOldToNewCompartmentId[comp->id] = n;
        comp->id = n;
        n++;
    }
    if (mechInstanceMap)
    {
        p.at(0) = 0; //top node gets parent Id 0 as in Coreneuron
        for (int i=1; i<p.size(); i++)
            p.at(i) = fromOldToNewCompartmentId.at(p.at(i));
    }

    ////// Mechanisms instances: merge all instances of all compartments into instances of the branch
    vector<vector<floble_t>> dataMechs (mechanismsCount);
    vector<vector<offset_t>> pdataMechs (mechanismsCount);
    vector<vector<offset_t>> nodesIndicesMechs (mechanismsCount);
    vector<vector<unsigned char>> vdataMechs (mechanismsCount);
    vector<deque<int>> mechsInstancesIds (mechanismsCount); //mech-offset -> list of pointers to mech instance value

    map< pair<int, offset_t>, offset_t> ionInstanceToDataOffset; //from pair of < ion mech type, OLD node id> to ion offset in NEW representation

    for (Compartment * comp : compartments)
    {
        int compDataOffset=0;
        int compPdataOffset=0;
        int compVdataOffset=0;
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
            if (mech->pntMap > 0 || mech->vdataSize>0)
            {
                assert(( type == IClamp && mech->vdataSize == 1 && mech->pdataSize == 2 && mech->pntMap>0)
                    || ( type == StochKv && mech->vdataSize == 1 && mech->pdataSize == 5 && mech->pntMap==0)
                    || ((type == ProbAMPANMDA_EMS || type == ProbGABAAB_EMS)
                        && mech->vdataSize == 2 && mech->pdataSize == 3 && mech->pntMap>0 ));

                size_t totalVdataSize = 0;
                if (type == IClamp || type == ProbAMPANMDA_EMS || type == ProbGABAAB_EMS)
                    totalVdataSize += sizeof(Point_process);
                if (type == StochKv || type == ProbAMPANMDA_EMS || type == ProbGABAAB_EMS)
                    totalVdataSize += sizeof(nrnran123_State);
                vdataMechs[mechOffset].insert(vdataMechs[mechOffset].end(), &comp->vdata[compVdataOffset], &comp->vdata[compVdataOffset+totalVdataSize] );
                compVdataOffset += totalVdataSize;
            }
        }
    }

    //merge all mechanisms vectors in the final one
    //store the offset of each mechanism data (for later)
    for (int m=0; m<mechanismsCount; m++)
    {
        Mechanism * mech = mechanisms[m];
        int dataOffset=0;
        int pdataOffset=0;
        int vdataOffset=0;
        for (int i=0; i<nodesIndicesMechs[m].size(); i++) //for all instances
        {
            assert(ionInstanceToDataOffset.find( make_pair(mech->type, nodesIndicesMechs[m][i]) ) == ionInstanceToDataOffset.end() );
            if (mech->isIon && inputParams->branchingDepth>0) //for pdata calculation
                ionInstanceToDataOffset[ make_pair(mech->type, nodesIndicesMechs[m][i]) ] = data.size();
            data.insert ( data.end(), &dataMechs[m][dataOffset], &dataMechs[m][dataOffset+mech->dataSize ]);
            nodesIndices.push_back(nodesIndicesMechs[m][i]);
            dataOffset  += mech->dataSize;

            if (mech->pntMap > 0 || mech->vdataSize>0)
            {
                int totalVdataSize=0;
                if (mech->type == IClamp || mech->type == ProbAMPANMDA_EMS || mech->type == ProbGABAAB_EMS)
                    totalVdataSize += sizeof(Point_process);
                if (mech->type == StochKv || mech->type == ProbAMPANMDA_EMS || mech->type == ProbGABAAB_EMS)
                    totalVdataSize += sizeof(nrnran123_State);
                assert(totalVdataSize>0);
                vdata.insert(vdata.end(), &vdataMechs[m][vdataOffset], &vdataMechs[m][vdataOffset+totalVdataSize]);
                vdataOffset += totalVdataSize;
            }

            if (inputParams->branchingDepth>0) //if we need to recalculate offsets or remove padding
            {
              for (int p=pdataOffset; p<pdataOffset+mech->pdataSize; p++)
              {
                offset_t pd = pdataMechs.at(m).at(p);
                int ptype = memb_func[mech->type].dparam_semantics[p-pdataOffset];
                switch (ptype)
                {
                case -1:  //"area" (6th field)
                {
                    assert(pd>=N*5 && pd<N*6);
                    offset_t oldId = pd-N*5;
                    offset_t newId = fromOldToNewCompartmentId.at(oldId);
                    pdataMechs.at(m).at(p) = n*5+newId;
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
                    pdataMechs.at(m).at(p) = (offset_t) vdataPointerOffset++;
                    break;
                case -8: //"bbcorepointer"
                    assert(0); //watch condition, not supported
                    break;
                default:
                    if (ptype>0 && ptype<1000) //name preffixed by '#'
                    {
                        //ptype is the ion (mechanism) type it depends on
                        //pdata is an offset in nt->data (a var in the ion)
                        //pd points to SoA notation, independently of the LAYOUT (converted before)

                        Mechanism * ion = neurox::GetMechanismFromType(ptype);
                        Mechanism::Ion ionOffset = ion->GetIonIndex();
                        IonInstancesInfo & ionInfo = ionsInstancesInfo.at((int) ionOffset);
                        int dataStart = ionInfo.dataStart;
                        assert(pd>=dataStart && pd<ionInfo.dataEnd);

                        int instanceOffset = trunc( (double)(pd-dataStart) / (double) ion->dataSize);
                        int instanceVariableOffset = (pd-dataStart) % ion->dataSize;
                        int nodeId = ionInfo.nodeIds.at(instanceOffset);
                        int newNodeId = fromOldToNewCompartmentId.at(nodeId);
                        pdataMechs.at(m).at(p) = ionInstanceToDataOffset.at(make_pair(ion->type, newNodeId)) + instanceVariableOffset;
                        assert(pdataMechs.at(m).at(p)>=n*6);
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
            }
            pdata.insert(pdata.end(), &pdataMechs[m][pdataOffset], &pdataMechs[m][pdataOffset+mech->pdataSize]);
            pdataOffset += mech->pdataSize;
        }
        dataMechs[m].clear();
        pdataMechs[m].clear();
        vdataMechs[m].clear();
        nodesIndicesMechs[m].clear();

        //convert neuron mech-instances ids from neuron- to branch-level
        if (mechInstanceMap)
          for (int i=0; i<mechsInstancesIds[m].size(); i++)
          {
            int oldInstanceId = mechsInstancesIds[m][i];
            (*mechInstanceMap)[m][oldInstanceId] = i;
          }
    }
    return n;
}

bool CompareCompartmentsPtrsIds(Compartment * a, Compartment *b)
{   return a->id<b->id; }

hpx_t DataLoader::CreateBranch(int nrnThreadId, hpx_t somaAddr, BranchType branchType, int thvar_index,
                               deque<Compartment*> & allCompartments, Compartment * topCompartment,
                               vector<DataLoader::IonInstancesInfo> & ionsInstancesInfo, int branchingDepth)
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
    int N = allCompartments.size();

    //Vector Play instances
    vector<floble_t> vecPlayT;
    vector<floble_t> vecPlayY;
    vector<PointProcInfo> vecPlayInfo;

    //branch NetCons
    vector<NetConX> branchNetCons;
    vector<neuron_id_t> branchNetConsPreId;
    vector<floble_t> branchWeights;

    hpx_t branchAddr = HPX_NULL; ///hpx address of this branch
    int thresholdVoffset=-1; //offset of the AP threshold voltage, or -1 if none in this branch
    Compartment * bottomCompartment = nullptr; //iterator for subsections

    if (inputParams->branchingDepth==0 ) //Flat a la Coreneuron
    {
        n = GetBranchData(allCompartments, data, pdata, vdata, p, instancesCount, nodesIndices, N, ionsInstancesInfo, NULL);
        GetVecPlayBranchData(allCompartments, vecPlayT, vecPlayY, vecPlayInfo, NULL);
        GetNetConsBranchData(allCompartments, branchNetCons, branchNetConsPreId, branchWeights, NULL);
        branchAddr = somaAddr;
        thresholdVoffset = thvar_index;
    }
    else if (inputParams->branchingDepth>0 ) //branch-parallelism
    {
        deque<Compartment*> subSection;
        if (branchingDepth>0)  //node of a tree, with branches
        {
            //subsection is the set of all sequential compartments until finding a bifurcation
            subSection.push_back(topCompartment);
            for (bottomCompartment = topCompartment;
                 bottomCompartment->branches.size()==1;
                 bottomCompartment = bottomCompartment->branches.front())
                subSection.push_back(bottomCompartment->branches.front());
        }
        else //leaf of the tree
        {
            //subsection is the set of all children compartments (recursively)
            GetAllChildrenCompartments(subSection, topCompartment);
        }

        //this step is only necesary so that data has the same alignment as CoreNeuron
        //and allows one to compare results . (can be removed for non-debug mode)
        std::sort(subSection.begin(), subSection.end(), CompareCompartmentsPtrsIds);

        //create sub-section of branch
        vector<map<int,int>> mechInstanceMap(mechanismsCount); //mech-offset -> ( map[old instance]->to new instance )
        GetMechInstanceMap(subSection, mechInstanceMap);
        n = GetBranchData(subSection, data, pdata, vdata, p, instancesCount, nodesIndices, N, ionsInstancesInfo, &mechInstanceMap);
        GetVecPlayBranchData(subSection, vecPlayT, vecPlayY, vecPlayInfo, &mechInstanceMap);
        GetNetConsBranchData(subSection, branchNetCons, branchNetConsPreId, branchWeights, &mechInstanceMap);

        //allocate children branches recursively (if any)
        if (bottomCompartment)
            for (size_t c=0; c<bottomCompartment->branches.size(); c++)
                branches.push_back(CreateBranch( nrnThreadId, somaAddr,
                    branchType==BranchType::Soma && c==0 ? BranchType::AxonInitSegment : BranchType::Dendrite,
                    branchType==BranchType::Soma && c==0 ? thvar_index - n : -1, /*offset in AIS = offset in soma - nt->end */
                    allCompartments, bottomCompartment->branches[c], ionsInstancesInfo, branchingDepth-1));

        //Benchmark and assign this branch to least busy compute node (except soma and AIS)
        //Note: we do this after children creation so that we use top (lighter) branches to balance work load
        hpx_t tempBranchAddr = hpx_gas_alloc_local(1, sizeof(Branch), NEUROX_MEM_ALIGNMENT);
        bool runBenchmarkAndClear = true;
        int dumbThresholdOffset=0;
        double timeElapsed=-1;
        hpx_call_sync(tempBranchAddr, Branch::Init,
                      &timeElapsed, sizeof(timeElapsed), //output
                      &n, sizeof(offset_t),
                      &nrnThreadId, sizeof(int),
                      &dumbThresholdOffset, sizeof(int),
                      data.size()>0 ? data.data() : nullptr, sizeof(floble_t)*data.size(),
                      pdata.size()>0 ? pdata.data() : nullptr, sizeof(offset_t)*pdata.size(),
                      instancesCount.data(), instancesCount.size()*sizeof(offset_t),
                      nodesIndices.data(), nodesIndices.size()*sizeof(offset_t),
                      &somaAddr, sizeof(hpx_t),
                      nullptr, 0, //no branching
                      p.data(), sizeof(offset_t)*p.size(), //force use of parent index
                      vecPlayT.size() > 0 ? vecPlayT.data() : nullptr, sizeof(floble_t)*vecPlayT.size(),
                      vecPlayY.size() > 0 ? vecPlayY.data() : nullptr, sizeof(floble_t)*vecPlayY.size(),
                      vecPlayInfo.size() > 0 ? vecPlayInfo.data() : nullptr, sizeof(PointProcInfo)*vecPlayInfo.size(),
                      branchNetCons.size() > 0 ? branchNetCons.data() : nullptr, sizeof(NetConX)*branchNetCons.size(),
                      branchNetConsPreId.size() > 0 ? branchNetConsPreId.data() : nullptr, sizeof(neuron_id_t)*branchNetConsPreId.size(),
                      branchWeights.size() > 0 ? branchWeights.data() : nullptr, sizeof(floble_t)*branchWeights.size(),
                      vdata.size()>0 ? vdata.data() : nullptr, sizeof(unsigned char)*vdata.size(),
                      &runBenchmarkAndClear, sizeof(bool)
                      );
        assert(timeElapsed>0);
        hpx_gas_clear_affinity(tempBranchAddr);

        //get HPX address of branch; create it if necessary, and update benchmark table
        int rank = hpx_get_my_rank();
        switch (branchType)
        {
          case BranchType::Soma :
            //already allocated
            branchAddr = somaAddr;
            //update benchmark table with this entry
            hpx_call_sync(HPX_THERE(0), tools::LoadBalancing::QueryLoadBalancingTable,
                          NULL, 0, //output
                          &timeElapsed, sizeof(double),
                          &rank, sizeof(int));
            break;
          case BranchType::AxonInitSegment :
            //AIS will be allocated in the same locality as soma (to speed-up AT threshold check)
            branchAddr = hpx_gas_alloc_local(1, sizeof(Branch), NEUROX_MEM_ALIGNMENT);
            thresholdVoffset = thvar_index; //correct value past by soma

            //update benchmark table with this entry
            hpx_call_sync(HPX_THERE(0), tools::LoadBalancing::QueryLoadBalancingTable,
                          NULL, 0, //output
                          &timeElapsed, sizeof(double),
                          &rank, sizeof(int));
            break;
          case BranchType::Dendrite :
            //get rank where to allocate it (query also updates benchmark table)
            hpx_call_sync(HPX_THERE(0), tools::LoadBalancing::QueryLoadBalancingTable,
                          &rank, sizeof(int), //output
                          &timeElapsed, sizeof(double));
            //allocate memory on that rank
            branchAddr = hpx_gas_alloc_local_at_sync(1, sizeof(Branch), NEUROX_MEM_ALIGNMENT, HPX_THERE(rank));

            //clean temporary memory
            break;
        }
#ifndef NDEBUG
        printf("-- %s branch of neuron nrn_id %d allocated to rank %d (%.6f ms)\n",
               branchType==BranchType::Soma ? "soma" : (branchType==BranchType::AxonInitSegment ? "AIS" : "dendrite"),
               nrnThreadId, rank, timeElapsed);
#endif
    }
    assert(branchAddr!= HPX_NULL);
    hpx_call_sync(branchAddr, Branch::Init,
                  NULL, 0, //no timing
                  &n, sizeof(offset_t),
                  &nrnThreadId, sizeof(int),
                  &thresholdVoffset, sizeof(int),
                  data.size()>0 ? data.data() : nullptr, sizeof(floble_t)*data.size(),
                  pdata.size()>0 ? pdata.data() : nullptr, sizeof(offset_t)*pdata.size(),
                  instancesCount.data(), instancesCount.size()*sizeof(offset_t),
                  nodesIndices.data(), nodesIndices.size()*sizeof(offset_t),
                  &somaAddr, sizeof(hpx_t),
                  branches.size() ? branches.data() : nullptr, branches.size() ? sizeof(hpx_t)*branches.size() : 0,
                  branches.size() ? nullptr : p.data(), branches.size() ? 0 : sizeof(offset_t)*p.size(),
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

hpx_action_t DataLoader::InitNetcons = 0;
int DataLoader::InitNetcons_handler()
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(DataLoader::InitNetcons);

    if (local->soma && inputParams->outputNetconsDot)
        fprintf(fileNetcons, "%d [style=filled, shape=ellipse];\n", local->soma->gid);

    std::deque<std::pair<hpx_t, floble_t> > netcons; //set of <srcAddr, minDelay> synapses to notify
    std::deque<std::pair<int, spike_time_t> > dependencies; //set of <srcGid, nextNotificationtime> for dependencies

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

            //tell the neuron to add a synapse to this branch and inform him of the fastest netcon we have
            if (minDelay != impossiblyLargeDelay) //if any active synapse
            {
                //add this netcon to list of netcons to be communicated
                netcons.push_back(make_pair(srcAddr,minDelay));

                //add this pre-syn neuron as my time-dependency
                if (inputParams->algorithm == AlgorithmType::All || inputParams->algorithm == AlgorithmType::BackwardEulerTimeDependencyLCO)
                {
                  spike_time_t notificationTime = inputParams->tstart+minDelay*
                          TimeDependencyLCOAlgorithm::TimeDependencies::notificationIntervalRatio;
                  dependencies.push_back( make_pair(srcGid, notificationTime ) );
                }
            }
        }
    }

    //inform pre-syn neuron that he connects to me
    hpx_t netconsLCO = hpx_lco_and_new(netcons.size());
    hpx_t topBranchAddr = local->soma ? target : local->branchTree->topBranchAddr;
    int myGid = local->soma ? local->soma->gid : -1;
    for (std::pair<hpx_t, floble_t> & nc : netcons)
        hpx_call(nc.first, DataLoader::AddSynapse, netconsLCO,
                 &target, sizeof(hpx_t), &nc.second, sizeof(nc.second),
                 &topBranchAddr, sizeof(hpx_t), &myGid, sizeof(int));
    hpx_lco_wait_reset(netconsLCO); hpx_lco_delete_sync(netconsLCO);

    //inform my soma of my time dependencies
    hpx_t dependenciesLCO = hpx_lco_and_new(dependencies.size());
    bool initPhase = true;
    for (std::pair<int, spike_time_t> & dep : dependencies)
        hpx_call(topBranchAddr, Branch::UpdateTimeDependency, dependenciesLCO,
                 &dep.first, sizeof(neuron_id_t), &dep.second, sizeof(spike_time_t),
                 &initPhase, sizeof(bool));
    hpx_lco_wait_reset(dependenciesLCO); hpx_lco_delete_sync(dependenciesLCO);

    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t DataLoader::AddSynapse = 0;
int DataLoader::AddSynapse_handler(
        const int nargs, const void *args[], const size_t[] )
{
    neurox_hpx_pin(Branch);
    assert(local->soma);
    assert(nargs==4);
    hpx_t addr = *(const hpx_t*) args[0];
    floble_t minDelay = *(const floble_t*) args[1];
    hpx_t topBranchAddr = *(const hpx_t*) args[2];
    int destinationGid = *(const int*) args[3];
    local->soma->AddSynapse(new Neuron::Synapse(addr,minDelay,topBranchAddr,destinationGid));
    neurox_hpx_unpin;
}

void DataLoader::RegisterHpxActions()
{
    neurox_hpx_register_action(neurox_zero_var_action,     DataLoader::Init);
    neurox_hpx_register_action(neurox_zero_var_action,     DataLoader::InitMechanisms);
    neurox_hpx_register_action(neurox_zero_var_action,     DataLoader::InitNeurons);
    neurox_hpx_register_action(neurox_zero_var_action,     DataLoader::InitNetcons);
    neurox_hpx_register_action(neurox_zero_var_action,     DataLoader::Finalize);
    neurox_hpx_register_action(neurox_several_vars_action, DataLoader::AddSynapse);
    neurox_hpx_register_action(neurox_several_vars_action, DataLoader::AddNeurons);
    neurox_hpx_register_action(neurox_several_vars_action, DataLoader::UpdateMechanismsDependencies);
}
