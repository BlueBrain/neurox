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

#include "neurox/neurox.h"

using namespace std;

void CoreNeuronDataLoader::loadCoreNeuronData(int argc, char ** argv)
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

void CoreNeuronDataLoader::loadData(int argc, char ** argv)
{ 
    loadCoreNeuronData(argc, argv); //CoreNeuron

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

    //create the tree structure for all neurons and mechanism
    vector<Mechanism> mechanisms;
    for_each(nrn_threads, nrn_threads+nrn_nthread, [&](NrnThread & nt) {

        vector<Compartment> compartments;
        //reconstructs tree with solver values
        for (int n=nt.end-1; n>=0; n--)
        {
            Compartment & compartment = compartments[n];
            compartment.setSolverValues(nt._actual_a[n], nt._actual_b[n], nt._actual_d[n], nt._actual_v[n], nt._actual_rhs[n], nt._actual_area[n]);

            if ( n>=nt.ncell) //if it is not top node, i.e. has a parent
            {
                Compartment & parentCompartment = compartments[nt._v_parent_index[n]];
                parentCompartment.addChild(&compartment);
            }
        }

        //reconstructs mechanisms for compartments
        int dataOffset=0, pdataOffset=0;
        //int synOffset=0;
        for (NrnThreadMembList* tml = nt.tml; tml!=nullptr; tml = tml->next)
        {
            int mechId = tml->index;
            //Data unique to each mechanism type
            mechanisms[mechId] = Mechanism(nrn_prop_param_size_[mechId],
                    nrn_prop_dparam_size_[mechId],tml->ndependencies,
                    pnt_map[mechId], nrn_is_artificial_[mechId], tml->dependencies);

            //Mechanisms' application to each compartment
            Memb_list *& ml = tml->ml; //list of compartments this mechanism belongs to
            int dataSize  = mechanisms[mechId].dataSize;
            int pdataSize = mechanisms[mechId].pdataSize;
            for (int n=0; n<ml->nodecount; n++)
            {
                Compartment & compartment = compartments[ml->nodeindices[n]];
                compartment.addMechanism(mechId, &ml->data[dataOffset], dataSize, &ml->pdata[pdataOffset], pdataSize);
                dataOffset  += dataSize;
                pdataOffset += pdataSize;
            }

            /*
            //Incoming Synapses
            if (mechanisms[mechId].pntMap > 0) //if point process
            {
                for (int n=0; n<ml->nodecount; n++)
                {
                    //Compartment & source = compartments[ml->nodeindices[n]];
                    Point_process * sourcePntProc = &nt.pntprocs[synOffset+n];
                    NetCon & nc = nt.netcons[synOffset+n];
                    Point_process * targetPntProc = nc.target_;
                    nc.
                    //pntprocs->;
                }
                synOffset += ml->nodecount;
            }
            ///  OR   ////
            //traverse (Output)preSyns (local spike exchange / event delivery)
            for (int s=0; s<nt.n_presyn; s++)
            {
                int neuronId = nt.presyns[s]->gid_;
                InputPreSyn * ips = gid2in[preNeuronId];
                int offsetIps = ips->nc_index_;
                for (int i = 0; i<ips->nc_cnt_; i++)
                {
                    NetCon* ncPost = netcon_in_presyn_order_[offsetIps+i];
                    //ncPost->target_->
                }
            }*/
        }

        //We will create all neurons in HPX memory
        for (int n=0; n<nt.end; n++)
        {
            vector<Synapse> synapses;
            int neuronId = nt.presyns[n].gid_;
            PreSyn * outSynapses = gid2out.at(neuronId);
            double APThreshold = outSynapses->threshold_;
            for (int c=0; c<outSynapses->nc_cnt_; c++) //for all axonal contact (outgoing Synapses)
            {
                NetCon * nc = netcon_in_presyn_order_[outSynapses->nc_index_+c];
                Point_process * target = nc->target_;
                //TODO from the target-> vars we get the mech type, and the instance of that type.
                //we will need the hpx address of it
                //TODO mod files requires point process structure, should be added
                synapses.push_back(Synapse(*nc->weight_, nc->delay_, HPX_NULL)); //TODO HPX_NULL?
            }
            createNeuron(neuronId, compartments.at(n), mechanisms, APThreshold, synapses);
        }
    });

    int neuronsCount =  std::accumulate(nrn_threads, nrn_threads+nrn_nthread, 0, [](int n, NrnThread & nt){return n+nt.ncell;});
    createBrain(neuronsCount, mechanisms);

    nrn_cleanup(); //remode CoreNeuron data structs
}

void CoreNeuronDataLoader::createBrain(int neuronsCount, vector<Mechanism> & mechanisms)
{
    //create Brain HPX structure
    hpx_t neuronsAddr = hpx_gas_calloc_blocked(neuronsCount, sizeof(Neuron), NEUROX_HPX_MEM_ALIGNMENT);
    assert(neuronsAddr != HPX_NULL);
    int mechanismsCount = mechanisms.size();

    //serialize mechanisms dependencies
    int totalDependenciesCount = accumulate(mechanisms.begin(), mechanisms.end(), 0, [](int n, Mechanism & m){return n+m.dependenciesCount;});
    int * mechsDependencies = new int [totalDependenciesCount];
    int mechsOffset=0;
    for_each(mechanisms.begin(), mechanisms.end(), [&] (Mechanism & m)
        { memcpy(&mechsDependencies[mechsOffset], m.dependencies, sizeof(int)*m.dependenciesCount);
          mechsOffset+= m.dependenciesCount;
        });

    //Broadcasts (serialized) brain
    printf("Broadcasting Circuit info...\n");
    int e = hpx_bcast_rsync(Brain::initialize, neuronsCount, neuronsAddr, mechanisms.data(), mechanisms.size(), mechsDependencies);
    assert(e == HPX_SUCCESS);

    delete [] mechsDependencies;
}

void CoreNeuronDataLoader::createNeuron(int gid, Compartment & topCompartment, vector<Mechanism> & mechanisms, double APThreshold, vector<Synapse> & synapses)
{
    hpx_t topBranch = createBranch(&topCompartment, mechanisms);
    hpx_call_sync(brain->getNeuronAddr(gid), Neuron::initialize, NULL, 0, gid, topBranch, APThreshold, synapses.data(), synapses.size());
}

hpx_t CoreNeuronDataLoader::createBranch(Compartment * topCompartment, vector<Mechanism> & mechanisms)
{
    Compartment *comp = topCompartment;
    Branch b;
    b.n =0; //number of compartments
    for (Compartment *comp = topCompartment;
         comp->children.size()==1;
         comp = comp->children.front())
    {
        b.n++;
        //TODO copy vars and mechs
    };

    b.childrenCount = comp->children.size();
    if (b.childrenCount > 0)
    {
        b.children = new hpx_t[b.childrenCount];
        for (int c=0; c<b.childrenCount; c++)
            b.children[c]=createBranch(comp->children[c], mechanisms);
    }
    else b.children = nullptr;

    //Allocate HPX Branch
    hpx_t branchAddr = hpx_gas_calloc_local(1, sizeof(Branch), NEUROX_HPX_MEM_ALIGNMENT);
    hpx_call_sync(branchAddr, Branch::initialize, NULL, 0, &b, sizeof(Branch));
    return branchAddr;
}
