#include <stdio.h>
#include <string>
#include <set>
#include <map>
#include <list>
#include <tuple>
#include <algorithm>

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
#include "neurox/datatypes/Mechanism.h"
#include "neurox/datatypes/Branch.h"

using namespace std;
using namespace neurox::Input;
using namespace neurox::Input::Coreneuron;

FILE *fileNetcons;
std::set<int> allNeuronsIdsSet;

int DataLoader::getNeuronIdFromNrnThreadId(int nrn_id)
{
    return nrn_threads[nrn_id].presyns[0].gid_;
}

void DataLoader::compareDataStructuresWithCoreNeuron(Branch * branch)
{
#ifndef NDEBUG
    if (!branch->soma) return; //run only once
    NrnThread & nt = nrn_threads[0];
    assert(nt.end == branch->n);
    for (int i=0; i<branch->n; i++)
    {
        assert(nt._actual_a[i] == branch->a[i]);
        assert(nt._actual_b[i] == branch->b[i]);
        assert(nt._actual_d[i] == branch->d[i]);
        assert(nt._actual_v[i] == branch->v[i]);
        assert(nt._actual_rhs[i] == branch->rhs[i]);
        assert(nt._actual_area[i] == branch->area[i]);
        if (branch->p)
        {  assert(nt._v_parent_index[i] == branch->p[i]); }
    }

    //make sure that morphology data is correctly aligned in mem
    for (int i=0; i<6*branch->n; i++)
            {   assert(nt._data[i]==branch->data[i]); }

    int mechCount=0;
    int vdataOffset=0;
    for (NrnThreadMembList* tml = nt.tml; tml!=NULL; tml = tml->next) //For every mechanism
    {
        int type = tml->index;
        int m = mechanismsMap[type];
        Memb_list * ml = tml->ml; //Mechanisms application to each compartment
        Branch::MechanismInstance & instance = branch->mechsInstances[m];
        assert(ml->nodecount == instance.count);
        //assert(ml->_nodecount_padded == instance.instancesCount);
        int dataSize = mechanisms[m]->dataSize;
        int pdataSize = mechanisms[m]->pdataSize;
        for (int n=0; n<ml->nodecount; n++) //for every mech instance
        {
            assert(ml->nodeindices[n]==instance.nodesIndices[n]);
            for (int i=0; i<dataSize; i++)
            {   assert(ml->data[i]==instance.data[i]); }

            for (int i=0; i<pdataSize; i++)
            {
                assert(ml->pdata[i] == instance.pdata[i]);
                assert(nt._data[ml->pdata[i]] == branch->data[instance.pdata[i]]);
            }

            /* We comment this because it runs for NULL presyn
            if (mechanisms[m]->pntMap)
            {
                //compare point_processes (index 1)
                Point_process * pp = (Point_process *) nt._vdata[vdataOffset+1];
                Point_process * pp2 = (Point_process *) branch->vdata[vdataOffset+1];
                assert(pp->_type == pp2->_type );
                assert(pp->_i_instance == pp2->_i_instance );
                assert(pp->_tid == pp2->_tid );
                assert(pp->_presyn == pp2->_presyn );
                vdataOffset+= mechanisms[m]->vdataSize;
            }
            */
        }
        mechCount++;
    }
    assert(mechCount==mechanismsCount);
#endif
}

void DataLoader::fromHpxToCoreneuronDataStructs(
        const void * branch_ptr, Memb_list & membList,
        NrnThread & nrnThread, int mechType)
{
    Branch * branch = (Branch*) branch_ptr;
    int m=mechanismsMap[mechType];
    Branch::MechanismInstance * mechsInstances = branch->mechsInstances;
    membList.data  = mechsInstances[m].data;
    membList.pdata = mechsInstances[m].pdata;
    membList.nodecount = mechsInstances[m].count;
    membList._nodecount_padded = membList.nodecount;
    membList.nodeindices = mechsInstances[m].nodesIndices;
    membList._thread = NULL; //TODO: ThreadDatum never used ?
    nrnThread.end = branch->n;
    nrnThread.ncell = 1;
    nrnThread.weights = NULL; //TODO: FOR NOW, until it crashes
    nrnThread.mapping = NULL; //TODO is it used?
    nrnThread._actual_a = branch->a;
    nrnThread._actual_b = branch->b;
    nrnThread._actual_d = branch->d;
    nrnThread._actual_rhs = branch->rhs;
    nrnThread._actual_v = branch->v;
    nrnThread._actual_area = branch->area;
    nrnThread._data = branch->data;
    //TODO is this worth is? to avoid race conditions, they use these shadow as intermediate values storing
    //Why not be a vector in vdata instead of Nt?
    //(i think is because in vdata you have pointer per instance, not per mech type)
    //does the RNG in vdata need to be one per instance or could be one per type?
    nrnThread._shadow_rhs = (mechType == ProbAMPANMDA_EMS || mechType==ProbGABAAB_EMS) ? new double[mechsInstances[m].count]() : NULL;
    nrnThread._shadow_d = (mechType==ProbAMPANMDA_EMS || mechType==ProbGABAAB_EMS) ? new double[mechsInstances[m].count]() : NULL;
    nrnThread._dt = dt;
    nrnThread._t = t;
    //TODO this field is only used inside the capacitance mod function in capac.c
    nrnThread.cj = inputParams->secondorder ?  2.0/inputParams->dt : 1.0/inputParams->dt;
    nrnThread._vdata = branch->vdata;
    //compareDataStructuresWithCoreNeuron(branch);
}

void DataLoader::coreNeuronInitialSetup(int argc, char ** argv)
{
    char prcellname[1024], filesdat_buf[1024];

    // initialise default coreneuron parameters
    //initnrn(); //part of GlobalInfo constructor

    // handles coreneuron configuration parameters
    cn_input_params input_params;

    // read command line parameters
    input_params.read_cb_opts( argc, argv );

    // set global variables for start time, timestep and temperature
    t = input_params.tstart;
    dt = input_params.dt;
    rev_dt = (int)(1./dt);
    celsius = input_params.celsius;

    // full path of files.dat file
    sd_ptr filesdat=input_params.get_filesdat_path(filesdat_buf,sizeof(filesdat_buf));

    // memory footprint after mpi initialisation
    //report_mem_usage( "After HPX_Init" );

    // reads mechanism information from bbcore_mech.dat
    mk_mech( input_params.datpath );

    //report_mem_usage( "After mk_mech" );

    // create net_cvode instance
    mk_netcvode();

    // One part done before call to nrn_setup. Other part after.
    if ( input_params.patternstim ) {
        nrn_set_extra_thread0_vdata();
    }

    //report_mem_usage( "Before nrn_setup" );

    //pass by flag so existing tests do not need a changed nrn_setup prototype.
    nrn_setup_multiple = input_params.multiple;
    nrn_setup_extracon = input_params.extracon;

    // reading *.dat files and setting up the data structures, setting mindelay
    nrn_setup( input_params, filesdat, nrn_need_byteswap );

    report_mem_usage( "After nrn_setup " );

    // Invoke PatternStim
    if ( input_params.patternstim ) {
        nrn_mkPatternStim( input_params.patternstim );
    }

    /// Setting the timeout
    nrn_set_timeout(200.);

    // show all configuration parameters for current run
    input_params.show_cb_opts();
}

void DataLoader::addNetConsForThisNeuron(int neuronId, int preNeuronId, int netconsCount, int netconsOffset, map< int, vector<NetConX*> > & netcons)
{
    for (int s = 0; s<netconsCount; s++)
    {
      NetCon* nc = netcon_in_presyn_order_[netconsOffset + s];
      int postNeuronId = getNeuronIdFromNrnThreadId(nc->target_->_tid);

      //if synapse is not for this neuron...
      if (postNeuronId!=neuronId) continue;

      int mechType = nc->target_->_type;
      netcons[preNeuronId].push_back(
          new NetConX(mechType, (int) nc->target_->_i_instance, nc->delay_,
                      nc->weight_, pnt_receive_size[mechType], nc->active_));
    }
}

void DataLoader::coreNeuronFakeSteps() //can be deleted
{
    //nrn_finitialize
    int i;
    NrnThread* _nt;
    dt2thread(-1.);
    nrn_thread_table_check();
    clear_event_queue();
    nrn_spike_exchange_init();
    nrn_play_init(); /* Vector.play */
        ///Play events should be executed before initializing events
    for (i=0; i < nrn_nthread; ++i) {
        nrn_deliver_events(nrn_threads + i); /* The play events at t=0 */
    }
          for (_nt = nrn_threads; _nt < nrn_threads + nrn_nthread; ++_nt) {
            for (i=0; i < _nt->end; ++i) {
            (_nt->_actual_v[(i)]) = inputParams->voltage;
            }
          }
    for (i=0; i < nrn_nthread; ++i) {
        nrn_ba(nrn_threads + i, BEFORE_INITIAL);
    }

    /* the memblist list in NrnThread is already so ordered */
    for (i=0; i < nrn_nthread; ++i) {
        NrnThread* nt = nrn_threads + i;
        NrnThreadMembList* tml;
        for (tml = nt->tml; tml; tml = tml->next) {
            mod_f_t s = memb_func[tml->index].initialize;
            if (s) {
                // TODO COMMENTED (*s)(nt, tml->ml, tml->index);
            }
        }
    }
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

hpx_action_t DataLoader::createNeuron = 0;
int DataLoader::createNeuron_handler(const int *i_ptr, const size_t)
{
        neurox_hpx_pin(Branch);
        int i=*i_ptr;

        //reconstructs neurons
        NrnThread & nt = nrn_threads[i];
        int neuronId = getNeuronIdFromNrnThreadId(i);

        //======= 1 - reconstructs matrix (solver) values =======

        deque<Compartment*> compartments;
        for (int n=0; n<nt.end; n++)
            compartments.push_back(
                new Compartment(n, nt._actual_a[n], nt._actual_b[n], nt._actual_d[n],
                                nt._actual_v[n], nt._actual_rhs[n], nt._actual_area[n],
                                nt._v_parent_index[n]));

        //reconstructs parents tree
        for (int n=1; n<nt.end; n++) //exclude top (no parent)
        {
            Compartment * parentCompartment = compartments.at(nt._v_parent_index[n]);
            parentCompartment->addChild(compartments.at(n));
        }

        if (inputParams->outputCompartmentsDot)
        {
          FILE *fileCompartments = fopen(string("compartment"+to_string(neuronId)+".dot").c_str(), "wt");
          fprintf(fileCompartments, "graph G_%d\n{ bgcolor=%s;  node [shape=cylinder];\n", neuronId, DOT_PNG_BACKGROUND_COLOR );
          printSubClustersToFile(fileCompartments, compartments.at(0)); //add subclusters
          for (auto c : compartments) //draw edges
            for (auto k : c->branches)
                fprintf(fileCompartments, "%d -- %d %s;\n", c->id, k->id, "");
          fprintf(fileCompartments, "}\n");
          fclose(fileCompartments);
        }

        //======= 2 - reconstructs mechanisms instances ========

        int vdataTotalOffset=0;
        int dataTotalOffset=nt.end*6;
        int pointProcTotalOffset=0;
        map<int, pair<int,int>> offsetToInstance; //map of data offset -> mech instance (mech Id, node Id)
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
                if (mech->isIon)
                    offsetToInstance[dataTotalOffset] = make_pair(type, ml->nodeindices[n]);
                double * data = &ml->data[dataOffset];
                int * pdata = &ml->pdata[pdataOffset];
                void ** vdata = &nt._vdata[vdataTotalOffset];
                Compartment * compartment = compartments.at(ml->nodeindices[n]);
                assert(compartment->id == ml->nodeindices[n]);

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
                    assert(pp->_presyn == NULL); //PS is always NULL? Maybe circuit too small --> no synapses
                    //presyns i think are only used by nrniv/netcvode.cpp mechfor net_event (not for our HPX use case)
                    //see netstim mod files (vdata[ppvar[1*STRIDE]] is used by artcell_net_send() ).
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
                dataTotalOffset   += mech->dataSize;
                dataOffset  += mech->dataSize;
                pdataOffset += mech->pdataSize;
                vdataTotalOffset += mech->vdataSize;
            }
        }

        //======= 3 - reconstruct NetCons =====================

        map< int, vector<NetConX*> > netcons ; //netcons per pre-synaptic neuron offset (not id)

        //get all incoming synapses from neurons in other MPI ranks (or NrnThread?)
        for (std::map<int, InputPreSyn*>::iterator syn_it = gid2in.begin();
             syn_it!=gid2in.end(); syn_it++)
        {
            int preNeuronId = syn_it->first;
            //we can also get it from netcon_srcgid[nrn_thread_it][netcon_it]
            InputPreSyn * ips = syn_it->second;
            addNetConsForThisNeuron(neuronId, preNeuronId, ips->nc_cnt_, ips->nc_index_, netcons);
        }

        //get all incoming synapses from neurons in the same MPI rank (or NrnThread?)
        for (std::map<int, PreSyn*>::iterator syn_it = gid2out.begin();
             syn_it!=gid2out.end(); syn_it++)
        {
            int preNeuronId = syn_it->first;
            PreSyn * ps = syn_it->second;
            addNetConsForThisNeuron(neuronId, preNeuronId, ps->nc_cnt_, ps->nc_index_, netcons);
        }

        if (inputParams->outputNetconsDot)
        {
          int netConsFromOthers=0;
          for (auto nc : netcons)
          {
            int gid = nc.first;
            if (allNeuronsIdsSet.find(gid) == allNeuronsIdsSet.end())
                netConsFromOthers++;
            else
            {
                double minDelay=99999;
                for (auto ncv : nc.second) //get minimum delay between neurons
                    minDelay = std::min(minDelay, ncv->delay);
                fprintf(fileNetcons, "%d -> %d [label=\"%d (%.2fms)\"];\n", gid , neuronId, nc.second.size(), minDelay);
            }
          }
          #if OUTPUT_NETCONS_DOT_FILE_INCLUDE_OTHERS==true
            if (netConsFromOthers>0)
              fprintf(fileNetcons, "%s -> %d [label=\"%d\" fontcolor=gray color=gray arrowhead=vee fontsize=12];\n", "others", neuronId, netConsFromOthers);
          #endif
        }

        //======= 4 - reconstruct VecPlayContinuous events =======
        for (int v=0; v<nt.n_vecplay; v++)
        {
            VecPlayContinuous *vec = (VecPlayContinuous*) nt._vecplay[v];
            //discover node, mechanism and data offset id that *pd points to
            double *pd = vec->pd_;
            PointProcInfo ppi;
            ppi.nodeId=-1;
            for (NrnThreadMembList* tml = nt.tml; tml!=NULL; tml = tml->next) //For every mechanism
            {
                int type = tml->index;
                Mechanism * mech = getMechanismFromType(type);
                Memb_list * ml = tml->ml;
                int dataOffset=0;
                for (int n=0; n<ml->nodecount; n++)
                {
                   // if is this mechanism and this instance
                   if (&ml->data[dataOffset] <= pd && pd < &ml->data[dataOffset+mech->dataSize])
                   {
                       ppi.nodeId = ml->nodeindices[n];
                       ppi.mechType = type;
                       ppi.mechInstance = n;
                       ppi.instanceDataOffset = pd - &ml->data[dataOffset];
                       ppi.size = vec->t_->size();
                       continue;
                   }
                   dataOffset += mech->dataSize;
                }
                if (ppi.nodeId!=-1) continue;
            }
            assert(ppi.nodeId != -1);
            compartments.at(ppi.nodeId)->addVecPlay(vec->t_->data(), vec->y_->data(), ppi);
        }

        //======= 5 - recursively create branches tree ===========

        double APthreshold = nrn_threads[i].presyns[0].threshold_;
        createBranch(target, compartments, compartments.at(0), netcons, (int) compartments.size(), offsetToInstance);
        hpx_call_sync(target, Branch::initSoma, NULL, 0,
                      &neuronId, sizeof(int), &APthreshold, sizeof(double));
        for (auto c : compartments)
            delete c;
        for (auto nc: netcons)
            for (auto nvc : nc.second)
                delete nvc;
    neurox_hpx_unpin;
}

void DataLoader::loadData(int argc, char ** argv)
{
    coreNeuronInitialSetup(argc, argv);

    //we will walk a bit with coreneuron
    coreNeuronFakeSteps();

    int neuronsCount = std::accumulate(nrn_threads, nrn_threads+nrn_nthread, 0, [](int n, NrnThread & nt){return n+nt.ncell;});

    allNeuronsIdsSet = std::set<int>();
    for (int i=0; i<neuronsCount; i++)
    {
        assert(nrn_threads[i].ncell == 1);
        allNeuronsIdsSet.insert(getNeuronIdFromNrnThreadId(i));
    }

    if (inputParams->outputCompartmentsDot)
    {
      for (int i=0; i<neuronsCount; i++)
      {
        int neuronId = getNeuronIdFromNrnThreadId(i);
        FILE *fileCompartments = fopen(string("compartment"+to_string(neuronId)+"_NrnThread.dot").c_str(), "wt");
        fprintf(fileCompartments, "graph G%d\n{  node [shape=cylinder];\n", neuronId );

        //for all nodes in this NrnThread
        NrnThread * nt = &nrn_threads[i];
        for (int n=nt->ncell; n<nt->end; n++)
            fprintf(fileCompartments, "%d -- %d;\n", nt->_v_parent_index[n], n);
        fprintf(fileCompartments, "}\n");
        fclose(fileCompartments);
      }
    }

    /** Reconstructs unique data related to each mechanism*
     * nargs=3 where:
     * args[0] = array of all mechanisms info
     * args[1] = array of all mechanisms depedencies (children on mechanisms execution tree)
     * args[2] = array of all mechanisms names (sym)
     */
    std::vector<Mechanism> mechsData;
    std::vector<int> mechsSuccessorsId;
    std::vector<char> mechsSym;

    //Different nrn_threads[i] have diff mechanisms sets; we will get the list of unique mechs types
    map<int, NrnThreadMembList*> uniqueMechs;
    for (int i=0; i<neuronsCount; i++)
        for (NrnThreadMembList* tml = nrn_threads[i].tml; tml!=NULL; tml = tml->next)
            if (uniqueMechs.find(tml->index)==uniqueMechs.end())
                uniqueMechs[tml->index] = tml;

    for (auto mech_it : uniqueMechs)
    {
        NrnThreadMembList *tml = mech_it.second;
        int type = tml->index;
        assert(type==mech_it.first);
        vector<int> successorsIds;
        int dependenciesCount;

        if (inputParams->multiMex)
        {
          std::set<int> dependenciesIds ; //set of 'unique' dependencies (ignores several dependencies between same mechs pair)
          for (int i=0; i<neuronsCount; i++)
            for (NrnThreadMembList* tml2 = nrn_threads[i].tml; tml2!=NULL; tml2 = tml2->next) //for every 2nd mechanism
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
          if ( tml->index == capacitance)  //capacitance is not part of graph
          {
            dependenciesCount=0;
          }
          else
          {
            //all except second element (the one after capacitance) have 1 dependency
            auto secondMech = uniqueMechs.begin(); //get first elem of map
            std::advance(secondMech, 1); //advance 1 position
            dependenciesCount = type==secondMech->first ? 0 : 1;
            //all except last one have a successor
            //auto successorMech = mech_it;
            auto successorMech = uniqueMechs.find(mech_it.first);
            std::advance(successorMech, 1); //the mech immediately after in the map
            if (successorMech != uniqueMechs.end())
                successorsIds.push_back(successorMech->first);
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

    if (inputParams->outputMechanismsDot)
    {
      FILE *fileMechs = fopen(string("mechanisms.dot").c_str(), "wt");
      fprintf(fileMechs, "digraph G\n{ bgcolor=%s;\n", DOT_PNG_BACKGROUND_COLOR);
      fprintf(fileMechs, "graph [ratio=0.3];\n", "start");
      fprintf(fileMechs, "%s [style=filled, shape=Mdiamond, fillcolor=beige];\n", "start");
      fprintf(fileMechs, "%s [style=filled, shape=Mdiamond, fillcolor=beige];\n", "end");
      fprintf(fileMechs, "\"%s (%d)\" [style=filled, fillcolor=beige];\n",
            getMechanismFromType(capacitance)->sym, capacitance);
      for (int m =0; m< mechanismsCount; m++)
      {
        Mechanism * mech = mechanisms[m];
        if (mech->pntMap > 0) //if is point process make it dotted
            fprintf(fileMechs, "\"%s (%d)\" [style=dotted];\n", mech->sym, mech->type);

        if (mech->dependenciesCount==0 && mech->type!=capacitance) //top mechanism
            fprintf(fileMechs, "%s -> \"%s (%d)\";\n", "start", mech->sym, mech->type);

        if (mech->successorsCount==0 && mech->type!= capacitance) //bottom mechanism
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
      fprintf(fileNetcons, "digraph G\n{ bgcolor=%s; layout=circo;\n", DOT_PNG_BACKGROUND_COLOR);
      //fprintf(fileNetcons, "digraph G\n{ bgcolor=%s;\n", DOT_PNG_BACKGROUND_COLOR);
      #if OUTPUT_NETCONS_DOT_FILE_INCLUDE_OTHERS==true
        fprintf(fileNetcons, "others [color=gray fontcolor=gray];\n");
      #endif
    }

    printf("neurox::createNeurons...\n");
    hpx_par_for_sync( [&] (int i, void*) -> int
    {  return hpx_call_sync(getNeuronAddr(i), DataLoader::createNeuron, NULL, 0, &i, sizeof(i));
    }, 0, neuronsCount, NULL);

    //all neurons have been created, every branch will inform pre-syn ids that they are connected
    assert(allNeuronsIdsSet.size() == neuronsCount);
    std::vector<int>   allNeuronsIdsVec  (neuronsCount);
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
                 allNeuronsIdsVec.data(),  allNeuronsIdsVec.size() *sizeof(int),
                 allNeuronsAddrVec.data(), allNeuronsAddrVec.size()*sizeof(hpx_t));
    hpx_lco_wait(lco_neurons);
    hpx_lco_delete(lco_neurons, NULL);

    if (inputParams->outputNetconsDot)
    {
      fprintf(fileNetcons, "}\n");
      fclose(fileNetcons);
    }
}

void DataLoader::getNetConsBranchData(
        deque<Compartment*> & compartments, map<int, vector<NetConX*> > & netcons,
        vector<NetConX> & branchNetCons, vector<int> & branchNetConsPreId,
        vector<double> & branchNetConsArgs)
{
    //get size and allocate it
    int ncCount=0, argsCount=0;
    for (auto nc : netcons)
        for (auto nc2 : nc.second)
        {
            ncCount++;
            argsCount += nc2->argsCount;
        }

    branchNetCons = vector<NetConX> (ncCount);
    branchNetConsPreId = vector<int> (ncCount);
    branchNetConsArgs = vector<double> (ncCount);

    int i=0;
    for (auto nc : netcons)
        for (auto nc2 : nc.second)
        {
           assert(nc2!=NULL);
           assert(nc2->args != NULL);
           memcpy(&branchNetCons.at(i), nc2, sizeof(NetConX)); //push_back will call destructor!
           memcpy(&branchNetConsPreId.at(i), &nc.first, sizeof(int));
           branchNetConsArgs.insert(branchNetConsArgs.end(), nc2->args, nc2->args + nc2->argsCount);
           branchNetCons[i].args=NULL; //needed or will be pointing to original netcons data and will call destructor
           i++;
        }
}

void DataLoader::getVecPlayBranchData(deque<Compartment*> & compartments, vector<double> & vecPlayTdata,
                                     vector<double> & vecPlayYdata, vector<PointProcInfo> & vecPlayInfo)
{
    for (auto comp : compartments)
    {
        vecPlayTdata.insert(vecPlayTdata.end(), comp->vecPlayTdata.begin(), comp->vecPlayTdata.end());
        vecPlayYdata.insert(vecPlayYdata.end(), comp->vecPlayYdata.begin(), comp->vecPlayYdata.end());
        vecPlayInfo.insert(vecPlayInfo.end(), comp->vecPlayInfo.begin(), comp->vecPlayInfo.end());
    }
}

int DataLoader::getBranchData(deque<Compartment*> & compartments, vector<double> & data, vector<int> & pdata, vector<void*> & vdata,
                                  vector<int> & p, vector<int> & instancesCount, vector<int> & nodesIndices, int totalN, map<int, pair<int,int>> & offsetToInstance)
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

    vector<vector<double>> dataMechs (mechanismsCount);
    vector<vector<int>> pdataMechs (mechanismsCount);
    vector<vector<int>> nodesIndicesMechs (mechanismsCount);
    vector<vector<void*>> vdataMechs (mechanismsCount);

    int n=0;
    vector<int> pdataType; //type of pdata offset per pdata entry
    map< pair<int,int>, int> instanceToOffset; //from pair of < mech type, OLD node id> to offset in NEW representation
    map<int,int> fromNeuronToBranchId; //map of branch id per compartment id
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
        for (int n=0; n<instancesCount[m]; n++)
          for (int d=0; d<mech->pdataSize; d++)
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
    for (int i=0; i<pdata.size(); i++)
    {
        int p = pdata.at(i);
        int ptype = pdataType.at(i);
        switch (ptype)
        {
        case 0: //not registered (its an *_ion type)
            //do nothing, its a flag (the 'iontype', see eion.c)
            break;
        case -1:  //"area" (6th field)
        {
            assert(p>=totalN*5 && p<totalN<6);
            int oldId = p-totalN*5;
            int newId = fromNeuronToBranchId[oldId];
            pdata[i] = n*5+newId;
            break;
        }
        case -2: //"iontype"
        case -3: //"cvodeieq"
        case -5: //"pointer"
            assert(0); //not used
            break;
        case -4: //"netsend"
        case -6: //"pntproc"
        case -7: //"bbcorepointer"
            pdata[i] = vdataOffset++;
            break;
        default:
            if (ptype>0 && ptype<1000) //name preffixed by '#'
            {   //ptype is the ion (mechanism) type it depends on
                //pdata[i] is an offset in pdata
                std::pair<int, int> mechNodePair;
                int oldOffset = -1;
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
                int newOffset = instanceToOffset.at(mechNodePair);
                assert(newOffset>=n*6);
                pdata[i] = newOffset + (pdata[i]-oldOffset); //'pdata[i]-oldOffset' is the offset on the data vector for that instance
                assert(pdata[i]>=n*6);
            }
            else if (ptype>=1000) //name not preffixed
            {
                pdata[i] = ptype-1000; //just a value (concentration) summed with 1000
            }
            else
                throw std::runtime_error("Unknown pdata type %d (FLAG3)\n");
            break;
        }
    }
    return n;
}

hpx_t DataLoader::createBranch(hpx_t target, deque<Compartment*> & compartments, Compartment * topCompartment,
                               map< int, vector<NetConX*> > & netcons, int totalN, map<int, pair<int,int>> & offsetToInstance)
{
    assert(topCompartment!=NULL);
    int n = -1; //number of compartments in branch
    vector<double> data; //compartments info (RHS, D, A, B, V, AREA)*n
    vector<int> pdata; //pointers to data
    vector<int> p; //parent nodes index
    vector<int> instancesCount (mechanismsCount);
    vector<int> nodesIndices;
    vector<hpx_t> branches;
    vector<void*> vdata; //TODO should go away at some point,

    //Vector Play instances
    vector<double> vecPlayT;
    vector<double> vecPlayY;
    vector<PointProcInfo> vecPlayInfo;

    //branch NetCons
    vector<NetConX> branchNetCons;
    vector<int> branchNetConsPreId;
    vector<double> branchNetConsArgs;

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
        for (int c=0; c<comp->branches.size(); c++)
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
                  &n, sizeof(int),
                  HPX_NULL, 0, //NOT USED
                  data.size()>0 ? data.data() : nullptr, sizeof(double)*data.size(),
                  pdata.size()>0 ? pdata.data() : nullptr, sizeof(int)*pdata.size(),
                  instancesCount.data(), instancesCount.size()*sizeof(int),
                  nodesIndices.data(), nodesIndices.size()*sizeof(int),
                  multiSplix ? branches.data() : nullptr, multiSplix ? sizeof(hpx_t)*branches.size() : 0,
                  multiSplix ? nullptr : p.data(), multiSplix ? 0 : sizeof(int)*p.size(),
                  vecPlayT.size() > 0 ? vecPlayT.data() : nullptr, sizeof(double)*vecPlayT.size(),
                  vecPlayY.size() > 0 ? vecPlayY.data() : nullptr, sizeof(double)*vecPlayY.size(),
                  vecPlayInfo.size() > 0 ? vecPlayInfo.data() : nullptr, sizeof(double)*vecPlayInfo.size(),
                  branchNetCons.size() > 0 ? branchNetCons.data() : nullptr, sizeof(NetCon)*branchNetCons.size(),
                  branchNetConsPreId.size() > 0 ? branchNetConsPreId.data() : nullptr, sizeof(int)*branchNetConsPreId.size(),
                  branchNetConsArgs.size() > 0 ? branchNetConsArgs.data() : nullptr, sizeof(double)*branchNetConsArgs.size(),
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

    const int   * neuronsIds  = (const int*)   args[0]; //list of all existing neurons ids
    const hpx_t * neuronsAddr = (const hpx_t*) args[1]; //list of all existing neurons addrs (same order as ids)

    //create index of unique pre-neurons ids
    std::set<int> netconsPreNeuronIds;
    for (auto nc : local->netcons)
        netconsPreNeuronIds.insert(nc.first);

    //inform pre-synaptic neurons (once) that we connect
    for (int i=0; i<sizes[0]/sizeof(int); i++) //for all existing neurons
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
