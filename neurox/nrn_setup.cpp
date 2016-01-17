#include <stdio.h>
#include <string>
#include <set>
#include <map>

#include "coreneuron/nrnconf.h"
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrnoc/nrnoc_decl.h"
#include "coreneuron/nrniv/netcon.h"
#include "coreneuron/nrniv/netcvode.h"
#include "coreneuron/nrniv/nrniv_decl.h"
#include "coreneuron/nrniv/nrn_assert.h"

#include "neurox/nrn_setup.h"
GlobalInfo * globalInfo;

using namespace std;

static hpx_action_t initGlobalInfo = 0;
static int initGlobalInfo_handler(const GlobalInfo *inputGlobalInfo, const size_t size)
{
    //Make sure message arrived correctly, and pin memory
    hpx_t target = hpx_thread_current_target();
    uint64_t *local = NULL;
    if (!hpx_gas_try_pin(target, (void**) &local))
        return HPX_RESEND;

    //do the work
    globalInfo = new GlobalInfo();
    memcpy(globalInfo, inputGlobalInfo, size);

    //unpin and return success
    hpx_gas_unpin(target);
    return HPX_SUCCESS;
}

static hpx_action_t initNeuronInfo = 0;
static int initNeuronInfo_handler(const Neuron * inputNeuron,  const size_t size)
{
    hpx_t neuron_addr = hpx_thread_current_target();
    Neuron * neuron = NULL;
    if (!hpx_gas_try_pin(neuron_addr, (void**) &neuron))
        return HPX_RESEND;

    //copy header information
    memcpy(neuron, inputNeuron, size);

    //allocate branches
    if (!globalInfo->multiSplit) //full neurons per CPU
        neuron->branches = hpx_gas_calloc_local  (neuron->branchesCount, sizeof (Branch), NEUROX_HPX_MEM_ALIGNMENT);
    else //branches spread across all nodes
        neuron->branches = hpx_gas_calloc_blocked(neuron->branchesCount, sizeof (Branch), NEUROX_HPX_MEM_ALIGNMENT);
    assert(neuron->branches != HPX_NULL);

    hpx_gas_unpin(neuron_addr);
    return HPX_SUCCESS;
}

//TODO: this is a hack: only 1 cpu reads some neruons and redistributes
void nrn_setup_hpx()
{
    //1. get total neurons count
    unsigned neuron_it=0;
    for (int k=0; k<nrn_nthread; k++)
    {
        FILE *dotfile = fopen(string("graph"+to_string(HPX_LOCALITY_ID)+"_"+to_string(k)+".dot").c_str(), "wt");
        fprintf(dotfile, "graph G%d\n{\n", k );

        //for all nodes in this NrnThread
        NrnThread * nt = &nrn_threads[k];
        for (int i=nt->ncell; i<nt->end; i++)
        {
            int p = nt->_v_parent_index[i];
            fprintf(dotfile, "%d -- %d;\n", p, i);
            if (nt->_v_parent_index[i] < nt->ncell) //TODO: is this a good soma check?
                neuron_it++;
        }
        fprintf(dotfile, "}\n");
        fclose(dotfile);
    }
/*
    GlobalInfo infoToBroadcast;
    infoToBroadcast.neuronsCount = neurons_it;
    infoToBroadcast.neuronsAddr  = hpx_gas_calloc_blocked(neurons_it, sizeof(Neuron), NEUROX_HPX_MEM_ALIGNMENT);
    infoToBroadcast.multiSplit   = HPX_LOCALITIES > neurons_it;

    //Broadcast input arguments
    printf("\nBroadcasting global info...\n");
    int e = hpx_bcast_rsync(initGlobalInfo, &infoToBroadcast, sizeof (GlobalInfo));
    assert(e == HPX_SUCCESS);

    //2. get branching structure (all maps take neuron hpx addr_t as key)
    map<hpx_t, deque<int>> parents; //given a neuron, return parents/P vector
    std::map< pair<hpx_t, int>, int> nodeToNeuron; //to which neuron a node belongs to

    neuron_it=0;
    for (int k=0; k<nrn_nthread; k++)
    {
        NrnThread * nt = &nrn_threads[k];
        for (int i=nt->ncell; i<nt->end; i++)
        {
            int p = nt->_v_parent_index[i];
            if (p < nt->cell) //if is a soma (no parent)
            {
                hpx_t neuron_addr = hpx_addr_add(globalInfo->neuronsAddr, sizeof(Neuron)*neuron_it, sizeof(Neuron));
                parent[neuron_addr][i]=-1;
                neuron_it++;
                nodeToNeuron[i] = neuron_addr;
            }
            else
            {
                hpx_t neuron_addr = nodeToNeuron[p]; //belongs to same neuron as parent
                parents[neuron_addr][i]=p;
            }
        }
    }

    neuron_it=0;
    for (int k=0; k<nrn_nthread; k++)
    {
        NrnThread * nt = &nrn_threads[k];

        std::map<hpx_t, int> neuronPerNode;
        map<hpx_t, map<int,int> > childrenPerNodePerNeuron;
        for (int i=nt->ncell; i<nt->end; i++)
        {
            int p = nt->_v_parent_index[i];

            //which neuron it belongs to?
            if ( p < nt->ncell) //if is soma
            {
                hpx_t neuron_addr = hpx_addr_add(globalInfo->neuronsAddr, sizeof(Neuron)*neuron_it, sizeof(Neuron));
                parentsIdsPerNeuron[neuron_addr].insert(p);
                nodesCountPerNeuron[neuron_addr]++;

                neuronPerNode[i] = neuron_addr;
            }
            else
            {
                hpx_t neuron_addr = neuronPerNode[p]; //same neuron as parent
                neuronPerNode[i] = neuron_addr;
                nodesCountPerNeuron[neuron_addr]++;

                //increment parent node's children count (if any)
                if (childrenPerNodePerNeuron[neuron_addr].find(p) == childrenPerNodePerNeuron[neuron_addr].end())
                    childrenPerNodePerNeuron[neuron_addr][p] = 1;
                else
                {
                    childrenPerNodePerNeuron[neuron_addr][p] ++;
                    parentsIdsPerNeuron[neuron_addr].insert(p);
                }
            }
        }
    }

    //allocate GAS memory for neuron and its branches
    printf("\nBroadcasting Neurons info...\n");
    hpx_addr_t done = hpx_lco_and_new(globalInfo->neuronsCount);
    for (int i=0; i<globalInfo->neuronsCount; i++)
    {
        hpx_t neuron_addr = hpx_addr_add(globalInfo->neuronsAddr, sizeof(Neuron)*t, sizeof(Neuron));
        Neuron neuron;
        neuron.branchesCount = parentsIdsPerNeuron[neuron_addr].size();
        neuron.neuronMetaData = 123;
        hpx_call(neuron_addr, initNeuronsInfo, done, &neuron, sizeof (Neuron));
    }
    hpx_lco_wait(done);
    hpx_lco_delete(done, HPX_NULL);

    for (int k=0; k<nrn_nthread; k++)
    {
        NrnThread * nt = &nrn_threads[k];
        for (int i=nt->ncell; i<nt->end; i++)
        {
            //TODO: what's the p of the very first node?
            int    p    = nt->_v_parent_index[i];
            double a    = nt->_actual_a[i];
            double b    = nt->_actual_b[i];
            double d    = nt->_actual_d[i];
            double rhs  = nt->_actual_rhs[i];
            double v    = nt->_actual_v[i];
            double area = nt->_actual_area[i];



            if (p<nt->ncell) //soma?
                somasIds.insert(p);
            if (p!= i-1) //
                children[p].push_back(i);
        }

*/
        /*
        double * nt_data = nt->_data; //nt->_ndata;
        //nidata, etc etc

        //for all mechanisms on this NrnThread
        for (NrnThreadMembList* tml = nt->tml; tml; tml = tml->next)
        {
            int type = tml->index;


            int is_art = nrn_is_artificial_[type];
            int layout = nrn_mech_data_layout_[type];
            char pntype = pnt_map[type];
            int is_ion = nrn_is_ion(type); //mem_func[type]...

            //Member functions (see declaration of initialize, etc)
            Memb_func* memfunc = &memb_func[type];

            ////TODO: dependencies not used on the code, right?
            //for (int d=0; d<tml->ndependencies; d++)
            //    dependenciesPerMech[tml].push_back( tml->dependencies[d]);

            Memb_list* ml = tml->ml;

            if (!is_art)
            for (int n=0; n< ml->nodecount; n++)
            {
                int nindex = ml->nodeindices[n];
                mechsPerNode[nindex].push_back(ml);
            }

            //data and its size
            int szp = nrn_prop_param_size_[type];
            double * data = ml->data;

            //pdata and its size
            int szdp = nrn_prop_dparam_size_[type];
            int * pdata = ml->pdata;

        }
        */

    /*
        //traverse preSyns (local spike exchange / event delivery)
        for (int s=0; s<nt->n_presyn; s++)
        {
            PreSyn * preSyn = &nt->presyns[s];
            for (int c=0; c<preSyn->nc_cnt_; c++)
            {
                NetCon * nc = netcon_in_presyn_order_[preSyn->nc_index_+c];
                Point_process * target = nc->target_;
                target->_tid; //TODO is this destination node id?
                            double * weight = nc->weight_;
                            bool active = nc->active_;
                            double delay = nc->delay_;
            }
        }

        //Remote spike exchange:
        //1. sender sends to a gid (cell id) i.e. a compute node
        //2. that compute node knows which node takes the incoming spike on that cell
        //(remote spike exchange / event delivery)

        //traverse all incoming events
        //gid2in returns the correspondence between sending gid and inputPreSyn
        for (std::map<int, InputPreSyn*>::iterator it = gid2in.begin(); it!= gid2in.end(); it++)
        {
            int sender = it->first;
            InputPreSyn * inputPreSyn = it->second;

            for (int c=0; c<inputPreSyn->nc_cnt_; c++)
            {
                NetCon * nc = netcon_in_presyn_order_[inputPreSyn->nc_index_+c];
                Point_process * target = nc->target_;
                target->_tid; //TODO is this destination node id?
                double * weight = nc->weight_;
                bool active = nc->active_;
                double delay = nc->delay_;
            }
        }

        nt->_ecell_memb_list; //??

        delete [] neurons;
    }

    */
}



void nrn_setup_register_hpx_actions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  initGlobalInfo, initGlobalInfo_handler, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  initNeuronInfo, initNeuronInfo_handler, HPX_POINTER, HPX_SIZE_T);
}

