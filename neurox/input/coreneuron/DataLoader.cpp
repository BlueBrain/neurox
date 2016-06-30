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

using namespace std;
using namespace Neurox;
using namespace Neurox::Input;

void CoreNeuronDataLoader::coreNeuronInitialSetup(int argc, char ** argv)
{
    //memory footprint after HPX initialisation
    report_mem_usage( "After hpx_init" );

    char prcellname[1024], filesdat_buf[1024], datpath[1024];

    // initialise default coreneuron parameters
    //initnrn(); //part of GlobalInfo constructor

    // handles coreneuron configuration parameters
    cn_input_params input_params;

    // read command line parameters
    input_params.read_cb_opts( argc, argv );

    // set global variables for start time, timestep and temperature
    t = input_params.tstart; // input_params.tstart;
    dt = input_params.dt ; //input_params.dt;
    rev_dt = 1/input_params.dt; //(int)(1./dt);
    celsius = input_params.celsius ; //input_params.celsius;

    // full path of files.dat file
    sd_ptr filesdat=input_params.get_filesdat_path(filesdat_buf,sizeof(filesdat_buf));

    // memory footprint after mpi initialisation
    report_mem_usage( "After MPI_Init" );

    // reads mechanism information from bbcore_mech.dat
    mk_mech( datpath );

    report_mem_usage( "After mk_mech" );

    // create net_cvode instance
    mk_netcvode();

    // One part done before call to nrn_setup. Other part after.
    if ( input_params.patternstim ) {
        nrn_set_extra_thread0_vdata();
    }

    report_mem_usage( "Before nrn_setup" );

    // reading *.dat files and setting up the data structures
    nrn_setup( input_params, filesdat, nrn_need_byteswap);

    report_mem_usage( "After nrn_setup " );

    // Invoke PatternStim
    if ( input_params.patternstim) {
        nrn_mkPatternStim( input_params.patternstim );
    }

    /// Setting the timeout
    nrn_set_timeout(200.);

    // show all configuration parameters for current run
    input_params.show_cb_opts();
}

//TODO this is a hack to access the ion_global_map var defined in eion.c only
//#include "coreneuron/nrnoc/eion.c"
//extern double** ion_global_map;
void CoreNeuronDataLoader::loadData(int argc, char ** argv)
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
    mechanismsCount = 0;
    for (NrnThreadMembList* tml = nrn_threads[0].tml; tml!=nullptr; tml = tml->next)
        mechanismsCount++;
    mechanisms = new Mechanism[mechanismsCount];

    for (NrnThreadMembList* tml = nrn_threads[0].tml; tml!=nullptr; tml = tml->next) //For every mechanism
    {
        int type = tml->index;
        mechanisms[type] = Mechanism(type, nrn_prop_param_size_[type],
                nrn_prop_dparam_size_[type],tml->ndependencies,
                pnt_map[type], nrn_is_artificial_[type], tml->dependencies, nrn_is_ion(type),
                -1, -1, -1);
                //ion_global_map[type][0], ion_global_map[type][1], ion_global_map[type][2]);
    }

    //requires one neuron per NRN_THREAD: in Blue config "CellGroupSize 1"
    neuronsCount =  std::accumulate(nrn_threads, nrn_threads+nrn_nthread, 0, [](int n, NrnThread & nt){return n+nt.ncell;});
    assert(neuronsCount == nrn_nthread);
    //allocate HPX memory space for neurons
    neuronsAddr = hpx_gas_calloc_cyclic(neuronsCount, sizeof(Neuron), NEUROX_HPX_MEM_ALIGNMENT);
    assert(neuronsAddr != HPX_NULL);

    //create map of incoming and outgoing synapses
    typedef struct SynapseOut
    {
        hpx_t synapseTargetAddr;
        int synapseTargetGid;
        double synapseDelay;
    };

    typedef struct SynapseIn
    {
        int preNeuronGid;
        int mechId;
        int mechInstance;
        double weight;
    };

    map<int, vector<SynapseOut> > outSynapses;
    map<int, vector<SynapseIn > > inSynapses;
    for (int i=0; i<nrn_nthread; i++)
    {
        NrnThread & nt = nrn_threads[i];
        int preNeuronId = nt.presyns[n].gid_;
        PreSyn * ps = gid2out[preNeuronId];
        //for all outgoing synapses...
        for (int i = 0; i<ps->nc_cnt; i++)
        {
            NetCon* nc = netcon_in_presyn_order_[ps->nc_index_ + i];

            int postNeuronId = nc->target_->_tid;

            SynapseIn synIn;
            synIn.weight = nc->weight_;
            synIn.mechType = nc->target_->_type;
            synIn.mechInstance = nc->target_->_i_instance;
            synIn.preNeuronGid = nc->target_->_tid;
            inSynapses[postNeuronId].push_back(synIn);

            SynapseOut synOut;
            synOut.synapseTargetAddr = HPX_NULL; //to be populated later
            synOut.synapseTargetGid = postNeuronId;
            synOut.synapseDelay = nc->delay_;
            outSynapses[preNeuronId].push_back(synOut);
        }
    }

    //create the tree structure for all neurons and mechanisms
    map<int,vector<Compartment>> compartments; //compartments per NrnThread Id
    map< tuple<int, int, int, int> , Compartment* > fromMechToCompartment;
    map< int, map <Compartment*, tuple<int, int, int, int> > > fromNeuronCompartmentToMech;

    for (int i=0; i<nrn_nthread; i++)
    {
        NrnThread & nt = nrn_threads[i];
        int neuronId = nt.presyns[n].gid_;

        //reconstructs tree with solver values
        for (int n=nt.end-1; n>=0; n--)
        {
            Compartment & compartment = compartments[neuronId][n];
            compartment.setSolverValues(nt._actual_a[n], nt._actual_b[n], nt._actual_d[n], nt._actual_v[n], nt._actual_rhs[n], nt._actual_area[n]);

            if ( n>=nt.ncell) //if it is not top node, i.e. has a parent
            {
                Compartment & parentCompartment = compartments[neuronId][nt._v_parent_index[n]];
                parentCompartment.addChild(&compartment);
            }
        }

        //reconstructs mechanisms for compartments
        int pdataOffset = 0;
        int dataOffset  = 6*nt.end; //a,b,d,v,rhs,area
        for (NrnThreadMembList* tml = nt.tml; tml!=nullptr; tml = tml->next) //For every mechanism
        {
            int type = tml->index;
            Memb_list *& ml = tml->ml; //Mechanisms application to each compartment
            int dataSize  = mechanisms[type].dataSize;
            int pdataSize = mechanisms[type].pdataSize;
            for (int n=0; n<ml->nodecount; n++) //for every compartment this mech type is applied to
            {
                Compartment & compartment = compartments[neuronId][ml->nodeindices[n]];
                compartment.addMechanism(type, n, &ml->data[dataOffset], dataSize, &ml->pdata[pdataOffset], pdataSize);
                dataOffset  += dataSize;
                pdataOffset += pdataSize;
            }
        }
        createNeuron(i, compartments[neuronId][0], gid2out[neuronId]->threshold_);
    }

    //reconstructs synapses
    //We will "spike" once every synapse and receive as return the hpx address of the branch as return value
    //First we build a map between synaptic info <NrnThread,type,n> to compartment
    for (int i=0; i<nrn_nthread; i++)
    {
        NrnThread & nt = nrn_threads[i];
        int neuronId = nt.presyns[n].gid_;
        for (NrnThreadMembList* tml = nt.tml; tml!=nullptr; tml = tml->next) //For every mechanism
        {
            int type = tml->index;
            Memb_list *& ml = tml->ml;
            for (int n=0; n<ml->nodecount; n++) //for every compartment this mech type is applied to
            {
                int compartmentId = ml->nodeindices[n];
                Compartment & compartment = compartments[neuronId][compartmentId];
                fromMechToCompartment[make_tuple(neuronId, type, n, compartmentId)] = &compartment;
            }
        }
    }

    for (int i=0; i<nrn_nthread; i++)
    {
        vector<hpx_t> outgoingSynapses;
        NrnThread & nt = nrn_threads[i];
        int neuronId = nt.presyns[n].gid_;
        if (gid2out.find(neuronId)==gid2out.end()) continue; //no outgoing synapses
        PreSyn * outSynapses = gid2out.at(preSynGid); //one PreSyn per neuron
        hpx_addr_t lco = outSynapses->nc_cnt_>0 ? hpx_lco_and_new(outSynapses->nc_cnt_) : HPX_NULL;
        for (int c=0; c<outSynapses->nc_cnt_; c++) //for all axonal contacts (outgoing Synapses)
        {
            NetCon * nc = netcon_in_presyn_order_[outSynapses->nc_index_+c];
            Point_process * target = nc->target_;
            int targetNrn = target->_tid; //or neuron if only 1 neuron per NrnThread
            int type = target->_type;
            int instance = target->_i_instance;

            tuple<int, int, int> keyTuple = make_tuple(targetNrn, type, instance);
            Compartment * comp = fromMechToCompartment.at(keyTuple);
            hpx_t synapseBranch = HPX_NULL;
            hpx_call_sync(getNeuronAddr(targetNrn), Neuron::addIncomingSynapse, &synapseBranch, sizeof(synapseBranch),
                          comp, *nc->weight_, nc->delay_, type, instance);
            outgoingSynapses.push_back(synapseBranch);
            /*
            //traverse tree up until find parent
            int postSynCompartmentId = nrn_threads[targetNrn]._ml_list[type]->nodeindices[instance];
            int postSynNeuronId = postSynCompartmentId;
            while (postSynNeuronId > nrn_threads[targetNrn].ncell)
                postSynNeuronId = nrn_threads[targetNrn]._v_parent_index[postSynNeuronId];
            */
        }
        if (lco != HPX_NULL)
        {
            hpx_lco_wait(lco);
            hpx_lco_delete(lco, HPX_NULL);

        }
        hpx_call_sync(getNeuronAddr(neuronId), Neuron::addOutgoingSynapses, nullptr, 0, outgoingSynapses.data(), outgoingSynapses.size());
    }
}

void CoreNeuronDataLoader::createNeuron(int gid, Compartment & topCompartment, double APthreshold)
{
    //recursively create tree
    hpx_t topBranch = createBranch(&topCompartment);
    //create neuron meta-data;
    hpx_call_sync(getNeuronAddr(gid), Neuron::init, NULL, 0, gid, topBranch, APthreshold);
}

hpx_t CoreNeuronDataLoader::createBranch(Compartment * topCompartment)
{
    int n=0; //number of compartments
    vector<double> d, b, a, rhs, v, area; //compartments info
    vector<int> mechsCount(mechanismsCount,0);
    vector<vector<double>> data(mechanismsCount); //data per mechanism id
    vector<vector<Datum>> pdata(mechanismsCount); //pdata per mechanism id
    Compartment *comp = nullptr;
    int m=0; //number of mechs instances

    for (comp = topCompartment;
         comp->branches.size()==1;
         comp = comp->branches.front())
    {
        d.push_back(comp->d);
        b.push_back(comp->b);
        a.push_back(comp->a);
        rhs.push_back(comp->rhs);
        v.push_back(comp->v);
        area.push_back(comp->area);

        //copy all mechanisms sorted by type
        for (int i=0; i<comp->mechsIds.size(); i++)
        {
            int mechId = comp->mechsIds[i];
            mechsCount[mechId]++;
            data[i].insert(data[i].end(), comp->data.begin(), comp->data.end());
            pdata[i].insert(pdata[i].end(), comp->pdata.begin(), comp->pdata.end());
            m++;
        }
        n++;
    };

    int branchesCount = comp->branches.size(); //children on final compartment of the branch
    hpx_t * branches = HPX_NULL;
    if (branchesCount > 0)
    {
        branches = new hpx_t[branchesCount];
        for (int c=0; c<branchesCount; c++)
            branches[c]=createBranch(comp->branches[c]);
    }

    //merge all data and pdata vectors into one
    vector<double> data_merged;
    vector<Datum> pdata_merged;
    data_merged.reserve(m);
    pdata_merged.reserve(m);
    for (int i=0; i<mechanismsCount; i++)
    {
        data_merged.insert(data_merged.end(), data[i].begin(), data[i].end() );
        pdata_merged.insert(pdata_merged.end(), pdata[i].begin(), pdata[i].end() );
    }

    //Allocate HPX Branch
    hpx_t branchAddr = hpx_gas_calloc_local(1, sizeof(Branch), NEUROX_HPX_MEM_ALIGNMENT);
    hpx_call_sync(branchAddr, Branch::init, NULL, 0, n, a.data(), b.data(), d.data(),
                  v.data(), rhs.data(), area.data(), m, mechsCount.data(), data_merged.data(),
                  pdata_merged.data(), branchesCount, branches);

    return branchAddr;
}
