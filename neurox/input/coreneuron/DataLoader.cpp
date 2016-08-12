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

#include "neurox/Neurox.h"
#include "neurox/datatypes/Mechanism.h"
#include "neurox/datatypes/Branch.h"

using namespace std;
using namespace NeuroX::Input;
using namespace NeuroX::Input::Coreneuron;

int DataLoader::getNeuronIdFromNrnThreadId(int nrn_id)
{
    nrn_threads[nrn_id].presyns[0].gid_;
}

void DataLoader::compareDataStructuresWithCoreNeuron(Branch * branch)
{
#ifdef DEBUG
    if (!branch->isSoma) return; //run only once
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
        void * branch_ptr, Memb_list & membList,
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
    nrnThread.cj = inputParams->secondorder ?  2.0/inputParams->dt : 1.0/inputParams->dt;
    //TODO shall this cj field be hardcoded on a branch info?
    nrnThread._vdata = branch->vdata;
    //Point_process ** pps = (Point_process **) nrn_threads[0]._vdata;

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
    report_mem_usage( "After HPX_Init" );

    // reads mechanism information from bbcore_mech.dat
    mk_mech( input_params.datpath );

    report_mem_usage( "After mk_mech" );

    // create net_cvode instance
    mk_netcvode();

    // One part done before call to nrn_setup. Other part after.
    if ( input_params.patternstim ) {
        nrn_set_extra_thread0_vdata();
    }

    report_mem_usage( "Before nrn_setup" );

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

void DataLoader::addNetConsForThisNeuron(int neuronId, int preNeuronId, int netconsCount, int netconsOffset, map<int, vector<NetConX*> > & netcons)
{
    for (int s = 0; s<netconsCount; s++)
    {
      NetCon* nc = netcon_in_presyn_order_[netconsOffset + s];
      int postNeuronId = getNeuronIdFromNrnThreadId(nc->target_->_tid);

      //if synapse is not for this neuron...
      if (postNeuronId!=neuronId) continue;

      int mechType = nc->target_->_type;
      NeuroX::NetConX * netcon = new NetConX(mechType, (int) nc->target_->_i_instance, nc->delay_,
                             nc->weight_, pnt_receive_size[mechType], nc->active_);
      netcons[preNeuronId].push_back(netcon);
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

void DataLoader::loadData(int argc, char ** argv)
{
    coreNeuronInitialSetup(argc, argv);

    //we will walk a bit with coreneuron
    coreNeuronFakeSteps();

#ifdef DEBUG
    for (int i=0; i<nrn_nthread; i++)
    {
        int neuronId = getNeuronIdFromNrnThreadId(i);

        FILE *fileCompartments = fopen(string("compartments"+to_string(neuronId)+"_NrnThread.dot").c_str(), "wt");
        fprintf(fileCompartments, "graph G%d\n{\n", i );

        //for all nodes in this NrnThread
        NrnThread * nt = &nrn_threads[i];
        for (int i=nt->ncell; i<nt->end; i++)
            fprintf(fileCompartments, "%d -- %d;\n", nt->_v_parent_index[i], i);
        fprintf(fileCompartments, "}\n");
        fclose(fileCompartments);
    }
#endif

    /** Reconstructs unique data related to each mechanism*
     * nargs=3 where:
     * args[0] = array of all mechanisms info
     * args[1] = array of all mechanisms depedencies (children on mechanisms execution tree)
     * args[2] = array of all mechanisms names (sym)
     */
    std::vector<Mechanism> mechsData;
    std::vector<int> mechsChildren;
    std::vector<char> mechsSym;
    for (NrnThreadMembList* tml = nrn_threads[0].tml; tml!=NULL; tml = tml->next) //For every mechanism
    {
        int type = tml->index;
        int symLength = memb_func[type].sym ? std::strlen(memb_func[type].sym) : 0;

        //TODO: for graph: remove circular or dual connections, (required mutex?)
        //no dependencies graph for now, next mechanism is the next in load sequence
        int childrenCount = tml->next==NULL ? 0 : 1;
        int children[1] = { tml->next==NULL ? -1 : tml->next->index };
        char isTopMechanism = type == capacitance ? 1 : 0; //exclude capacitance

        mechsData.push_back(
            Mechanism (type, nrn_prop_param_size_[type], nrn_prop_dparam_size_[type],
                       nrn_is_artificial_[type], pnt_map[type], nrn_is_ion(type),
                       symLength, NULL, //sym will be serialized below
                       isTopMechanism, childrenCount, NULL));  //children will be serialized below

        mechsChildren.insert(mechsChildren.end(), children, children + childrenCount);
        mechsSym.insert(mechsSym.end(), memb_func[type].sym, memb_func[type].sym + symLength);
    }

    printf("Broadcasting %d mechanisms...\n", mechsData.size());
    int e = hpx_bcast_rsync(NeuroX::setMechanisms,
                            mechsData.data(), sizeof(Mechanism)*mechsData.size(),
                            mechsChildren.data(), sizeof(int)* mechsChildren.size(),
                            mechsSym.data(), sizeof(char)*mechsSym.size());
    assert(e == HPX_SUCCESS);

    mechsData.clear(); mechsChildren.clear(); mechsSym.clear();

#ifdef DEBUG
    FILE *fileMechs = fopen(string("mechanisms.dot").c_str(), "wt");
    fprintf(fileMechs, "digraph G\n{\n");
    for (int m =0; m< mechanismsCount; m++)
    {
        Mechanism * mech = mechanisms[m];
        if (mech->isTopMechanism)
            fprintf(fileMechs, "%s -> %d;\n", "start", mech->type);

        for (int d=0; d<mech->childrenCount; d++)
            fprintf(fileMechs, "%d -> %d;\n", mech->type, mech->children[d]);
    }
    fprintf(fileMechs, "}\n");
    fclose(fileMechs);
#endif

    //requires one neuron per NRN_THREAD: in Blue config we must add: "CellGroupSize 1"
    int neuronsCount = std::accumulate(nrn_threads, nrn_threads+nrn_nthread, 0, [](int n, NrnThread & nt){return n+nt.ncell;});
    if(neuronsCount == nrn_nthread)
        printf("Warning: neurons count %d not equal to nrn_nthread %d\n", neuronsCount, nrn_nthread);

    //allocate HPX memory space for neurons
    printf("Broadcasting %d neurons...\n", neuronsCount);
    hpx_t neuronsAddr = hpx_gas_calloc_cyclic(neuronsCount, sizeof(Neuron), NEUROX_HPX_MEM_ALIGNMENT);
    e = hpx_bcast_rsync(NeuroX::setNeurons, &neuronsCount, sizeof(int), &neuronsAddr, sizeof(hpx_t));
    assert(e == HPX_SUCCESS);
    assert(neuronsAddr != HPX_NULL);

#ifdef DEBUG
    FILE *fileNetcons = fopen(string("netcons.dot").c_str(), "wt");
    std::set<int> availableNeuronsIds;
    for (int i=0; i<nrn_nthread; i++)
        availableNeuronsIds.insert(getNeuronIdFromNrnThreadId(i));
    fprintf(fileNetcons, "digraph G\n{\n");
#endif

    //reconstructs neurons
    for (int i=0; i<nrn_nthread; i++)
    {
        NrnThread & nt = nrn_threads[i];
        if (nt.ncell!=1)
        {
            printf("Warning: ignoring NrnThread %d because it has %d neurons instead of 1\n", i, nt.ncell);
            continue;
        }

        assert(nt.ncell==1); //only 1 cell per NrnThreadId;
        int neuronId =getNeuronIdFromNrnThreadId(i);

        //reconstructs matrix (solver) values
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

#ifdef DEBUG
        FILE *fileCompartments = fopen(string("compartments"+to_string(neuronId)+"_HPX.dot").c_str(), "wt");
        fprintf(fileCompartments, "graph G%d\n{\n", neuronId );
        for (auto c : compartments)
        {
            if (c->branches.size()>1)
                fprintf(fileCompartments, "%d [style=filled, fillcolor=beige];\n",  c->id);
            for (auto k : c->branches)
                fprintf(fileCompartments, "%d -- %d;\n", c->id, k->id);
        }
        fprintf(fileCompartments, "}\n");
        fclose(fileCompartments);
#endif

        //reconstructs mechanisms for compartments
        int vdataOffset=0;
        int ntPointProcOffset=0;
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
                double * data = &ml->data[dataOffset];
                int * pdata = &ml->pdata[pdataOffset];
                void ** vdata = &nt._vdata[vdataOffset];
                Compartment * compartment = compartments.at(ml->nodeindices[n]);
                assert(compartment->id == ml->nodeindices[n]);
#ifdef DEBUG
                if (mech->pntMap > 0) //debugging only
                {
                    assert((type==IClamp && mech->pdataSize==2)
                        ||((type==ProbAMPANMDA_EMS || type==ProbGABAAB_EMS) && mech->pdataSize==3));

                    //Offsets in pdata: 0: data, 1: point proc in nt._vdata, [2: rng in nt._vdata]
                    Point_process * pp = &nt.pntprocs[ntPointProcOffset++];
                    assert(nt._vdata[pdata[1]] == pp);
                    assert(pp->_presyn == NULL); //PS is always NULL? Maybe circuit too small --> no synapses
                    //presyns i think are only used by nrniv/netcvode.cpp mechfor net_event (not for our HPX use case)

                    if (mech->pdataSize > 2)
                    {
                        assert(pdata[1]+1 ==pdata[2]); //TODO no need to store offsets 1 and 2 if they are sequential
                        void * RNG = nt._vdata[pdata[2]];
                    }
                    for (int v=0; v<mech->vdataSize; v++)
                      compartment->vdata.push_back(vdata[v]);
                }
                else
                {
                    assert(mech->vdataSize==0);
                }
#endif
                compartment->addMechanismInstance(type, data, mech->dataSize, pdata,  mech->pdataSize);
                dataOffset  += mech->dataSize;
                pdataOffset += mech->pdataSize;
                vdataOffset += mech->vdataSize;
            }
        }

        map<int, vector<NetConX*> > netcons; //netcons per pre-synaptic neuron id

        //get all incoming synapses from neurons in other MPI ranks (or NrnThread?)
        for (std::map<int, InputPreSyn*>::iterator syn_it = gid2in.begin();
             syn_it!=gid2in.end(); syn_it++)
        {
            int preNeuronId = syn_it->first;
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

#ifdef DEBUG
        int netConsFromOthers=0;
        for (auto nc : netcons)
        {
            if (availableNeuronsIds.find(nc.first) == availableNeuronsIds.end())
                netConsFromOthers++;
            else
                fprintf(fileNetcons, "%d -> %d [label=\"%d\"];\n", nc.first, neuronId, nc.second.size());
        }
        if (netConsFromOthers>0)
                fprintf(fileNetcons, "%s -> %d [label=\"%d\"];\n", "others", neuronId, netConsFromOthers);
#endif

        //recursively create morphological tree, neuron metadata, and synapses
        double APthreshold = nrn_threads[i].presyns[0].threshold_;
        //double APthreshold = gid2out.at(neuronId)->threshold_;
        hpx_t topBranch = createBranch((char) 1, compartments, compartments.at(0), netcons);
        hpx_call_sync(getNeuronAddr(neuronId), Neuron::init, NULL, 0,
                      &neuronId, sizeof(int),
                      &topBranch, sizeof(hpx_t),
                      &APthreshold, sizeof(double));

        for (auto c : compartments)
            delete c;

        for (auto nc_it : netcons)
            for (auto nc : nc_it.second)
                delete nc;
    }

#ifdef DEBUG
    fprintf(fileNetcons, "}\n");
    fclose(fileNetcons);
#endif
}

int DataLoader::getBranchData(deque<Compartment*> & compartments, vector<double> & data, vector<int> & pdata, vector<void*> & vdata,
                                  vector<int> & p, vector<int> & instancesCount, vector<int> & nodesIndices)
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
    for (auto comp : compartments)
    {
        int compDataOffset=0;
        int compPdataOffset=0;
        int compVdataOffset=0;
        for (int m=0; m<comp->mechsTypes.size(); m++) //for all instances
        {
            int mechType = comp->mechsTypes[m];
            int mechOffset = mechanismsMap[mechType];
            assert(mechOffset>=0 && mechOffset<mechanismsCount);
            Mechanism * mech = mechanisms[mechOffset];
            dataMechs[mechOffset].insert(dataMechs[mechOffset].end(), &comp->data[compDataOffset], &comp->data[compDataOffset+mech->dataSize] );
            pdataMechs[mechOffset].insert(pdataMechs[mechOffset].end(), &comp->pdata[compPdataOffset], &comp->pdata[compPdataOffset+mech->pdataSize] );
            vdataMechs[mechOffset].insert(vdataMechs[mechOffset].end(), &comp->vdata[compVdataOffset], &comp->vdata[compVdataOffset+mech->vdataSize] );
            nodesIndicesMechs[mechOffset].push_back(n);
            instancesCount[mechOffset]++;
            compDataOffset += mech->dataSize;
            compPdataOffset += mech->pdataSize;
            compVdataOffset += mech->vdataSize;
        }
        n++;
    }

    //merge all mechanisms vectors in the final one
    for (int m=0; m<mechanismsCount; m++)
    {
        data.insert(data.end(), dataMechs[m].begin(), dataMechs[m].end());
        pdata.insert(pdata.end(), pdataMechs[m].begin(), pdataMechs[m].end());
        vdata.insert(vdata.end(), vdataMechs[m].begin(), vdataMechs[m].end());
        nodesIndices.insert(nodesIndices.end(), nodesIndicesMechs[m].begin(), nodesIndicesMechs[m].end());
        dataMechs[m].clear();
        pdataMechs[m].clear();
        vdataMechs[m].clear();
        nodesIndicesMechs[m].clear();
    }
    return n;
}

hpx_t DataLoader::createBranch(char isSoma, deque<Compartment*> & compartments, Compartment * topCompartment,  map<int, vector<NetConX*> > & netcons)
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

    if (multiSpliX)
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
        n = getBranchData(subSection, data, pdata, vdata, p, instancesCount, nodesIndices);

        //recursively create children branches
        for (int c=0; c<comp->branches.size(); c++)
            branches.push_back(createBranch((char) 0, compartments, comp->branches[c], netcons));
    }
    else //Flat a la Coreneuron
    {
        n = getBranchData(compartments, data, pdata, vdata, p, instancesCount, nodesIndices);
    }

    //Allocate HPX Branch
    hpx_t branchAddr = hpx_gas_calloc_local(1, sizeof(Branch), NEUROX_HPX_MEM_ALIGNMENT);

    //initiate main branch information (n, a, b, d, v, rhs, area, branching and mechanisms)
    hpx_call_sync(branchAddr, Branch::init, NULL, 0,
                  &n, sizeof(int),
                  &isSoma, sizeof(char),
                  data.size()>0 ? data.data() : nullptr, sizeof(double)*data.size(),
                  pdata.size()>0 ? pdata.data() : nullptr, sizeof(int)*pdata.size(),
                  instancesCount.data(), instancesCount.size()*sizeof(int),
                  nodesIndices.data(), nodesIndices.size()*sizeof(int),
                  multiSpliX ? branches.data() : nullptr, multiSpliX ? sizeof(hpx_t)*branches.size() : 0,
                  multiSpliX ? nullptr : p.data(), multiSpliX ? 0 : sizeof(int)*p.size(),
                  vdata.size()>0 ? vdata.data() : nullptr, sizeof(void*)*vdata.size());

    return branchAddr;
}
