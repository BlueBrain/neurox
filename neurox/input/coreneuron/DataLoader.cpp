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

void DataLoader::fromHpxToCoreneuronDataStructs(
        Branch * branch, Memb_list & membList,
        NrnThread & nrnThread, int mechType)
{
    Branch::MechanismInstance * mechsInstances = branch->mechsInstances;
    membList.data  = mechsInstances[mechType].data;
    membList.pdata = mechsInstances[mechType].pdata;
    membList.nodecount = mechsInstances[mechType].instancesCount;
    membList.nodeindices = mechsInstances[mechType].nodesIndices;
    membList._thread = NULL; //TODO: ThreadDatum never used ?
    nrnThread._actual_d = branch->d;
    nrnThread._actual_rhs = branch->rhs;
    nrnThread._actual_v = branch->v;
    nrnThread._actual_area = branch->area;
    nrnThread._dt = dt;
    nrnThread._t = t;
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

void DataLoader::loadData(int argc, char ** argv)
{
    coreNeuronInitialSetup(argc, argv);

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

        //TODO: for graph: remove circular or dual connections, should be a regular tree not a graph
        //no dependencies graph for now, next mechanism is the next in load sequence
        int childrenCount = tml->next!=NULL ? 1 : 0;
        int children[1] = { tml->next ? (short) tml->next->index: (short) -1};
        char isTopMechanism = tml == nrn_threads[0].tml ? 1 : 0;

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
    int neuronsCount =  std::accumulate(nrn_threads, nrn_threads+nrn_nthread, 0, [](int n, NrnThread & nt){return n+nt.ncell;});
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
        vector<Compartment*> compartments;
        for (int n=0; n<nt.end; n++)
            compartments.push_back(
                new Compartment(n, nt._actual_a[n], nt._actual_b[n], nt._actual_d[n],
                                nt._actual_v[n], nt._actual_rhs[n], nt._actual_area[n],
                                nt._v_parent_index[n]));

        //reconstructs parents tree
        for (int n=nt.end-1; n>0; n--) //exclude top (no parent)
        {
            Compartment * parentCompartment = compartments[nt._v_parent_index[n]];
            parentCompartment->addChild(compartments[n]);
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
        for (NrnThreadMembList* tml = nt.tml; tml!=NULL; tml = tml->next) //For every mechanism
        {
            int pdataOffset = 0;
            int dataOffset  = 0; //6*nt.end; //(a,b,d,v,rhs,area for NrnThread->_data)
            int type = tml->index;
            Memb_list * ml = tml->ml; //Mechanisms application to each compartment
            Mechanism * mech = getMechanismFromType(type);
            int dataSize  = mech->dataSize;
            int pdataSize = mech->pdataSize;
            for (int n=0; n<ml->nodecount; n++) //for every mech instance (or compartment this mech is applied to)
            {
                //TODO: I think ml->data is vectorized!
                assert (ml->nodeindices[n] < compartments.size());
                Compartment * compartment = compartments[ml->nodeindices[n]];
                compartment->addMechanismInstance(type, n, &ml->data[dataOffset], dataSize, &ml->pdata[pdataOffset], pdataSize);
                dataOffset  += dataSize;
                pdataOffset += pdataSize;
            }
        }

        map<int, vector<NetConX*> > netcons; //netcons per pre-synaptic neuron id

        //get all incoming synapses from neurons in other MPI ranks
        for (std::map<int, InputPreSyn*>::iterator syn_it = gid2in.begin();
             syn_it!=gid2in.end(); syn_it++)
        {
            int preNeuronId = syn_it->first;
            InputPreSyn * ips = syn_it->second;
            addNetConsForThisNeuron(neuronId, preNeuronId, ips->nc_cnt_, ips->nc_index_, netcons);
        }

        //get all incoming synapses from neurons in the same MPI rank
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
        hpx_t topBranch = createBranch((char) 1, compartments.at(0), netcons);
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

Compartment * DataLoader::getBranchSectionData(Compartment * topCompartment, int & n, vector<double> & d, vector<double> & b,
                          vector<double> & a, vector<double> & rhs, vector<double> & v, vector<double> & area,
                          vector<int> & p, vector<int> & instancesCount, vector<vector<double>> & data,
                          vector<vector<int>> & pdata, vector<vector<int>> & nodesIndices, char multiSplit)
{
    //iterate through all compartments on the branch
    Compartment *comp = topCompartment;
    while(1)
    {
        d.push_back(comp->d);
        b.push_back(comp->b);
        a.push_back(comp->a);
        v.push_back(comp->v);
        rhs.push_back(comp->rhs);
        area.push_back(comp->area);
        p.push_back(comp->p);

        //copy all mechanisms instances
        int dataOffset=0;
        int pdataOffset=0;
        for (int m=0; m<comp->mechsTypes.size(); m++) //for all instances
        {
            int mechType = comp->mechsTypes[m];
            int mechOffset = mechanismsMap[mechType];
            assert(mechOffset>=0 && mechOffset<mechanismsCount);
            Mechanism * mech = mechanisms[mechOffset];

            data[mechOffset].insert(data[mechOffset].end(), &comp->data[dataOffset], &comp->data[dataOffset+mech->dataSize] );
            pdata[mechOffset].insert(pdata[mechOffset].end(), &comp->pdata[pdataOffset], &comp->pdata[pdataOffset+mech->pdataSize] );
            nodesIndices[mechOffset].push_back(n);
            dataOffset += mech->dataSize;
            pdataOffset += mech->pdataSize;
            instancesCount[mechOffset]++;
        }
        n++;

        if (comp->branches.size() == 0) //leaf
            return comp;

        if (comp->branches.size() > 1) //bifurcation
        {
            if (!multiSplit)
            {
              Compartment * bottomComp = NULL;
              for (int c=0; c<comp->branches.size(); c++)
                bottomComp = getBranchSectionData(comp->branches[c], n, d, b, a, rhs, v,
                                        area, p, instancesCount, data, pdata, nodesIndices,
                                        multiSplit);
            }
            return comp;
        }

        //otherwise, iterate (take the next compartment in the sequence)
        comp = comp->branches.front();
    }
    throw new std::runtime_error("Error while reconstructing morphology (FLAG1)\n");
}

hpx_t DataLoader::createBranch(char isSoma, Compartment * topCompartment,  map<int, vector<NetConX*> > & netcons)
{
    int n=0; //number of compartments
    vector<double> d, b, a, rhs, v, area; //compartments info
    vector<int> p; //parent nodes index

    vector<int> instancesCount (mechanismsCount); //instances count per mechanism type (initialized to 0)
    vector<vector<double>> data(mechanismsCount); //data per mechanism type
    vector<vector<int>> pdata(mechanismsCount); //pdata per mechanism type
    vector<vector<int>> nodesIndices(mechanismsCount); //nodes indices per mechanism type

    char multiSplit=0;
    Compartment * comp = getBranchSectionData(topCompartment, n, d, b, a, rhs, v,
                                              area, p, instancesCount, data, pdata, nodesIndices,
                                              multiSplit);

    vector<hpx_t> branches (comp->branches.size());

    if (multiSplit) //next branches will be *hpx children* of this one
    {
      //recursively create children branches
      for (int c=0; c<comp->branches.size(); c++)
        branches[c]=createBranch((char) 0, comp->branches[c], netcons);
    }

    //Allocate HPX Branch
    hpx_t branchAddr = hpx_gas_calloc_local(1, sizeof(Branch), NEUROX_HPX_MEM_ALIGNMENT);

    //initiate main branch information (n, a, b, d, v, rhs, area, branchesCount, branches)
    hpx_call_sync(branchAddr, Branch::init, NULL, 0,
                  &isSoma, sizeof(char),
                  a.data(), sizeof(double)*a.size(),
                  b.data(), sizeof(double)*b.size(),
                  d.data(), sizeof(double)*d.size(),
                  v.data(), sizeof(double)*v.size(),
                  rhs.data(), sizeof(double)*rhs.size(),
                  area.data(), sizeof(double)*area.size(),
                  multiSplit ? branches.data() : NULL, multiSplit ? sizeof(hpx_t)*branches.size() : 0,
                  multiSplit ? NULL : p.data(), multiSplit ? 0 : sizeof(int)*p.size());

    //merge all mechanisms vectors into the first one
    for (int m=1; m<mechanismsCount; m++)
    {
        data[0].insert(data[0].end(), data[m].begin(), data[m].end());
        pdata[0].insert(pdata[0].end(), pdata[m].begin(), pdata[m].end());
        nodesIndices[0].insert(nodesIndices[0].end(), nodesIndices[m].begin(), nodesIndices[m].end());
        data[m].clear();
        pdata[m].clear();
        nodesIndices[m].clear();
    }

    hpx_call_sync(branchAddr, Branch::initMechanismsInstances, NULL, 0,
                  instancesCount.data(), instancesCount.size()*sizeof(int),
                  data[0].data(), data[0].size()*sizeof(double),
                  pdata[0].data(), pdata[0].size()*sizeof(int),
                  nodesIndices[0].data(), nodesIndices[0].size()*sizeof(int));

    //serialize all netcons and initialize them
    //... TODO get mech instance righ,t and neurons id (preNeuronId is not in the circuit)
    //hpx_call_sync(branchAddr, Branch::initMechanismsInstances, NULL, 0);

    return branchAddr;
}
