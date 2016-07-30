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
using namespace Neurox::Input;
using namespace Neurox::Input::Coreneuron;

int DataLoader::getNeuronIdFromNrnThreadId(int nrn_id)
{
    nrn_threads[nrn_id].presyns[0].gid_;
}

void DataLoader::fromHpxToCoreneuronDataStructs(
        Branch * branch, Memb_list & membList,
        NrnThread & nrnThread, short int mechType)
{
    Branch::MechanismInstances * mechs = branch->mechsInstances;
    membList.data  = mechs[mechType].data;
    membList.pdata = mechs[mechType].pdata;
    membList.nodecount = mechs[mechType].instancesCount;
    membList.nodeindices = mechs[mechType].nodesIndices;
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

void DataLoader::addNetConsForThisNeuron(int neuronId, int preNeuronId, int netconsCount, int netconsOffset, map<int, vector<NetConX> > & netcons)
{
    for (int s = 0; s<netconsCount; s++)
    {
      NetCon* nc = netcon_in_presyn_order_[netconsOffset + s];
      int postNeuronId = getNeuronIdFromNrnThreadId(nc->target_->_tid);

      //if synapse is not for this neuron...
      if (postNeuronId!=neuronId) continue;

      short int mechType = (short int) nc->target_->_type;
      Neurox::NetConX netcon(mechType, (int) nc->target_->_i_instance, nc->delay_,
                             nc->weight_, pnt_receive_size[mechType], nc->active_);
      netcons[preNeuronId].push_back(netcon);
      netcons[preNeuronId].shrink_to_fit();
    }
}


//TODO this is a hack to access the ion_global_map var defined in eion.c only
//#include "coreneuron/nrnoc/eion.c"
//extern double** ion_global_map;
void DataLoader::loadData(int argc, char ** argv)
{ 
    coreNeuronInitialSetup(argc, argv);

    //TODO: Debug: plot morphologies as dot file
    for (int k=0; k<nrn_nthread; k++)
    {
        FILE *dotfile = fopen(string("graph"+to_string(HPX_LOCALITY_ID)+"_"+to_string(k)+".dot").c_str(), "wt");
        fprintf(dotfile, "graph G%d\n{\n", k );

        //for all nodes in this NrnThread
        NrnThread * nt = &nrn_threads[k];
        for (int i=nt->ncell; i<nt->end; i++)
            fprintf(dotfile, "%d -- %d;\n", nt->_v_parent_index[i], i);
        fprintf(dotfile, "}\n");
        fclose(dotfile);
    }  

    //Reconstructs unique data related to each mechanism
    /**
     * nargs=3 where:
     * args[0] = array of all mechanisms info
     * args[1] = array of all mechanisms dependencies
     * args[2] = array of all mechanisms names (sym)
     */
    std::vector<Mechanism> mechsData;
    std::vector<int> mechsDependencies;
    std::vector<char> mechsSym;
    for (NrnThreadMembList* tml = nrn_threads[0].tml; tml!=NULL; tml = tml->next) //For every mechanism
    {
        int type = tml->index;
        int symLength = memb_func[type].sym ? std::strlen(memb_func[type].sym) : 0;
        mechsData.push_back(
            Mechanism (type, nrn_prop_param_size_[type], nrn_prop_dparam_size_[type],
                       nrn_is_artificial_[type], pnt_map[type], nrn_is_ion(type),
                       symLength, NULL,
                       tml->ndependencies, NULL));

        mechsDependencies.insert(mechsDependencies.end(), tml->dependencies, tml->dependencies + tml->ndependencies);
        mechsSym.insert(mechsSym.end(), memb_func[type].sym, memb_func[type].sym + symLength);
    }

    printf("Broadcasting %d mechanisms...\n", mechsData.size());
    int e = hpx_bcast_rsync(Neurox::setMechanisms,
                            mechsData.data(), sizeof(Mechanism)*mechsData.size(),
                            mechsDependencies.data(), sizeof(int)* mechsDependencies.size(),
                            mechsSym.data(), sizeof(char)*mechsSym.size());
    assert(e == HPX_SUCCESS);

    //requires one neuron per NRN_THREAD: in Blue config we must add: "CellGroupSize 1"
    neuronsCount =  std::accumulate(nrn_threads, nrn_threads+nrn_nthread, 0, [](int n, NrnThread & nt){return n+nt.ncell;});
    assert(neuronsCount == nrn_nthread);

    //allocate HPX memory space for neurons
    neuronsAddr = hpx_gas_calloc_cyclic(neuronsCount, sizeof(Neuron), NEUROX_HPX_MEM_ALIGNMENT);
    assert(neuronsAddr != HPX_NULL);

    //reconstructs neurons
    for (int i=0; i<nrn_nthread; i++)
    {
        NrnThread & nt = nrn_threads[i];
        assert(nt.ncell==1); //only 1 cell per NrnThreadId;
        int neuronId =getNeuronIdFromNrnThreadId(i);

        //reconstructs compartments tree with solver values
        vector<Compartment> compartments; //compartments
        for (int n=nt.end-1; n>=0; n--)
        {
            Compartment & compartment = compartments[n];
            compartment.setValues(n, nt._actual_a[n], nt._actual_b[n], nt._actual_d[n], nt._actual_v[n], nt._actual_rhs[n], nt._actual_area[n]);

            if ( n>=nt.ncell) //if it is not top node, i.e. has a parent
            {
                Compartment & parentCompartment = compartments[nt._v_parent_index[n]];
                parentCompartment.addChild(&compartment);
            }
        }

        //reconstructs mechanisms for compartments
        int pdataOffset = 0;
        int dataOffset  = 6*nt.end; //a,b,d,v,rhs,area
        for (NrnThreadMembList* tml = nt.tml; tml!=NULL; tml = tml->next) //For every mechanism
        {
            int type = tml->index;
            Memb_list *& ml = tml->ml; //Mechanisms application to each compartment
            int dataSize  = mechanisms[type].dataSize;
            int pdataSize = mechanisms[type].pdataSize;
            for (int n=0; n<ml->nodecount; n++) //for every compartment this mech type is applied to
            {
                Compartment & compartment = compartments[ml->nodeindices[n]];
                compartment.addMechanism(type, n, &ml->data[dataOffset], dataSize, &ml->pdata[pdataOffset], pdataSize);
                dataOffset  += dataSize;
                pdataOffset += pdataSize;
            }
        }

        map<int, vector<NetConX>> netcons; //netcons per pre-synaptic neuron id

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

        //recursively create morphological tree, neuron metadata, and synapses
        double APthreshold = gid2out.at(neuronId)->threshold_;
        hpx_t topBranch = createBranch(&compartments.at(0), netcons);
        hpx_call_sync(getNeuronAddr(neuronId), Neuron::init, NULL, 0, neuronId, topBranch, APthreshold);
    }
}

hpx_t DataLoader::createBranch(Compartment * topCompartment,  map<int, vector<NetConX> > & netcons)
{
    int n=0; //number of compartments
    vector<double> d, b, a, rhs, v, area; //compartments info
    vector<int> instancesCount(mechanismsCount,0); //instances count per mechanism id
    vector<vector<double>> data(mechanismsCount); //data per mechanism id
    vector<vector<Datum>> pdata(mechanismsCount); //pdata per mechanism id
    vector<vector<int> > nodesIndices(mechanismsCount); //nodes indices per mechanisms id
    vector<hpx_t> branches; //children branches

    Compartment *comp = NULL;
    for (comp = topCompartment;
         comp->branches.size()==1;
         comp = comp->branches.front())
    {
        d.push_back(comp->d);
        b.push_back(comp->b);
        a.push_back(comp->a);
        v.push_back(comp->v);
        rhs.push_back(comp->rhs);
        area.push_back(comp->area);

        //copy all mechanisms instances (sorted by type)
        for (int m=0; m<comp->mechsIds.size(); m++)
        {
            int mechId = comp->mechsIds[m];
            instancesCount[mechId]++;
            data[mechId].insert (data[mechId].end(),  comp->data.begin(),  comp->data.end() );
            pdata[mechId].insert(pdata[mechId].end(), comp->pdata.begin(), comp->pdata.end());
            nodesIndices[mechId].push_back(n);
        }
        n++;
    }

    for (int c=0; c<comp->branches.size(); c++)
        branches[c]=createBranch(comp->branches[c], netcons);

    //Allocate HPX Branch (TODO single node implementation)
    hpx_t branchAddr = hpx_gas_calloc_local(1, sizeof(Branch), NEUROX_HPX_MEM_ALIGNMENT);
    hpx_call_sync(branchAddr, Branch::init, NULL, 0,
                  n, &a, &b, &d, &v, &rhs, &area, &instancesCount, &data,
                  &pdata, &nodesIndices, &branches, &netcons);

    return branchAddr;
}
