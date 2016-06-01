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

extern double** ion_global_map;
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

    //create the tree structure for all neurons and mechanism
    vector<Mechanism> mechanisms; //unique info for mechanism
    vector<vector<Compartment>> compartments; //compartments per NrnThread Id
    vector<vector<tuple<int,Synapse>>> synapses; //incoming synapses per neuron
        //for a post syn neuron, list of <preId, Synapse>

    //from nrn id, mech type, mech instance, to compartment
    map< tuple<int, int, int> , Compartment *> fromMechToCompartment;
    for (int i=0; <nrn_nthread; i++)
    {
        NrnThread & nt = nrn_threads[i];

        //reconstructs tree with solver values
        for (int n=nt.end-1; n>=0; n--)
        {
            Compartment & compartment = compartments[i][n];
            compartment.setSolverValues(nt._actual_a[n], nt._actual_b[n], nt._actual_d[n], nt._actual_v[n], nt._actual_rhs[n], nt._actual_area[n]);

            if ( n>=nt.ncell) //if it is not top node, i.e. has a parent
            {
                Compartment & parentCompartment = compartments[i][nt._v_parent_index[n]];
                parentCompartment.addChild(&compartment);
            }
        }

        //reconstructs mechanisms for compartments
        int dataOffset=0, pdataOffset=0;
        for (NrnThreadMembList* tml = nt.tml; tml!=nullptr; tml = tml->next) //For every mechanism
        {
            int type = tml->index;
            char isIon = nrn_is_ion(type);

            //Data unique to each mechanism type
            mechanisms[type] = Mechanism(type, nrn_prop_param_size_[type],
                    nrn_prop_dparam_size_[type],tml->ndependencies,
                    pnt_map[type], nrn_is_artificial_[type], tml->dependencies,
                    isIon, ion_global_map[type][0], ion_global_map[type][1], ion_global_map[type][2]);

            //Mechanisms application to each compartment
            Memb_list *& ml = tml->ml; //list of compartments this mechanism is applied to to
            int dataSize  = mechanisms[type].dataSize;
            int pdataSize = mechanisms[type].pdataSize;
            for (int n=0; n<ml->nodecount; n++) //for every compartment this mech type is applied to
            {
                Compartment & compartment = compartments[ml->nodeindices[n]];
                compartment.addMechanism(type, n, &ml->data[dataOffset], dataSize, &ml->pdata[pdataOffset], pdataSize);
                dataOffset  += dataSize;
                pdataOffset += pdataSize;
                //will allow us to reconstruct synapses
                fromMechToCompartment[make_tuple(i, type, n)] = &compartment;
            }
        }
    }

    //reconstructs synapses
    for (int i=0; i<nrn_nthread; i++)
    {
        NrnThread & nt = nrn_threads[i];
        for (int n=0; n<nt.ncell; n++)
        {
            int preSynGid = nt.presyns[n].gid_;
            PreSyn * outSynapses = gid2out.at(preSynGid); //one PreSyn per neuron

            if (outSynapses==gid2out.end()) continue; //no outcoming synapses

            for (int c=0; c<outSynapses->nc_cnt_; c++) //for all axonal contacts (outgoing Synapses)
            {
                NetCon * nc = netcon_in_presyn_order_[outSynapses->nc_index_+c];
                Point_process * target = nc->target_;
                int targetNrn = target->_tid;
                int type = target->_type;
                int instance = target->_i_instance;

                Compartment * comp = fromMechToCompartment(make_tuple(targetNrn, type, instance));
                comp->addSynapse(*nc->weight_, nc->delay_, type, instance);

                //traverse tree up until find parent
                int postSynCompartmentId = nrn_threads[targetNrn]._ml_list[type]->nodeindices[instance];
                int postSynNeuronId = postSynCompartmentId;
                while (postSynNeuronId > nrn_threads[targetNrn].ncell)
                    postSynNeuronId = nrn_threads[targetNrn]._v_parent_index[postSynGid];

                //Synapse synapse(postSynCompartmentId, *nc->weight_, nc->delay_, mechType, mechInstance);
                //synapses[postSynGid].push_back(std::make_tuple(preSynGid, synapse));
            }
        }
    }

    //Build neurons
    for (int i=0; i<nrn_nthread; i++)
    {
        NrnThread & nt = nrn_threads[i];
        for (int n=0; n<nt.ncell; n++)
        {
            int neuronId = nt.presyns[n].gid_;
            double APThreshold = outSynapses->threshold_;
            createNeuron(neuronId, compartments.at(i).at(n), mechanisms, APThreshold, synapses);
        }
    }

    //TODO do a collective call to convert Synapse addresses
    //from ThreadId, mechType, mechInstance (in NrnThread)
    //to branch hpx_t, mechType, mechInstance (in branch)

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
    int e = hpx_bcast_rsync(Brain::init, neuronsCount, neuronsAddr, mechanisms.data(), mechanisms.size(), mechsDependencies);
    assert(e == HPX_SUCCESS);

    delete [] mechsDependencies;
}

void CoreNeuronDataLoader::createNeuron(int gid, Compartment & topCompartment, vector<Mechanism> & mechanisms, double APThreshold, vector<Synapse> & synapses)
{
    hpx_t topBranch = createBranch(&topCompartment, mechanisms);
    hpx_call_sync(getNeuronAddr(gid), Neuron::init, NULL, 0, gid, topBranch, APThreshold, synapses.data(), synapses.size());
}

hpx_t CoreNeuronDataLoader::createBranch(Compartment * topCompartment, vector<Mechanism> & mechanisms)
{
    int mechsTypesCount = mechanisms.size();
    int n=0; //number of compartments
    vector<double> d, b, a, rhs, v, area; //compartments info
    vector<int> mechsCount(mechsTypesCount,0);
    vector<vector<double>> data(mechsTypesCount); //data per mechanism id
    vector<vector<Datum>> pdata(mechsTypesCount); //pdata per mechanism id
    Compartment *comp = nullptr;
    int m=0; //number of mechs
    for (comp = topCompartment;
         comp->children.size()==1;
         comp = comp->children.front())
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

    int childrenCount = comp->children.size(); //final compartment of the branch
    hpx_t * children = HPX_NULL;
    if (childrenCount > 0)
    {
        children = new hpx_t[childrenCount];
        for (int c=0; c<childrenCount; c++)
            children[c]=createBranch(comp->children[c], mechanisms);
    }

    //merge all vectors into one
    vector<double> data_merged;
    vector<Datum> pdata_merged;
    data_merged.reserve(m);
    pdata_merged.reserve(m);
    for (int i=0; i<mechsTypesCount; i++)
    {
        data_merged.insert(data_merged.end(), data[i].begin(), data[i].end() );
        pdata_merged.insert(pdata_merged.end(), pdata[i].begin(), pdata[i].end() );
    }

    //Allocate HPX Branch
    hpx_t branchAddr = hpx_gas_calloc_local(1, sizeof(Branch), NEUROX_HPX_MEM_ALIGNMENT);
    hpx_call_sync(branchAddr, Branch::init, NULL, 0, n, a.data(), b.data(), d.data(),
                  v.data(), rhs.data(), area.data(), m, mechsCount.data(), data_merged.data(),
                  pdata_merged.data(), childrenCount, children);

    return branchAddr;
}
