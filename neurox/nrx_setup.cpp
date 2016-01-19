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

static hpx_action_t initializeGlobalInfo = 0;
static int initializeGlobalInfo_handler(const GlobalInfo *inputGlobalInfo, const size_t size)
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

static hpx_action_t initializeNeuron = 0;
static int initializeNeuron_handler(const Neuron * inputNeuron,  const size_t size)
{
    //Make sure message arrived correctly, and pin memory
    hpx_t neuron_addr = hpx_thread_current_target();
    Neuron * neuron = NULL;
    if (!hpx_gas_try_pin(neuron_addr, (void**) &neuron))
        return HPX_RESEND;

    //copy header information
    memcpy(neuron, inputNeuron, size);

    //unpin and return success
    hpx_gas_unpin(neuron_addr);
    return HPX_SUCCESS;
}

static hpx_action_t initializeBranch = 0;
static int initializeBranch_handler(const byte * branch_serial_input,  const size_t size)
{
    //Make sure message arrived correctly, and pin memory
    hpx_t branch_addr = hpx_thread_current_target();
    Branch * branch = NULL;
    if (!hpx_gas_try_pin(branch_addr, (void**) &branch))
        return HPX_RESEND;

    branch->deserialize(branch_serial_input, size);

    //unpin and return success
    hpx_gas_unpin(branch_addr);
    return HPX_SUCCESS;
}

hpx_t nrx_setup_branch(NrnThread * nt, map<int, list<int> > & tree, map<int, list< pair<int, Memb_list*> > > & mechanisms, int topNodeId )
{
    vector<double> a,b,d,rhs,v,area; //nodes information
    map<int,int> fromNodeToIndex; //map from general nodes ids to sequential ids

    int n=0, i=topNodeId;
    while (1)
    {
        fromNodeToIndex[i]=n++;

        //node level info
        a.push_back(nt->_actual_a[i]);
        b.push_back(nt->_actual_b[i]);
        d.push_back(nt->_actual_d[i]);
        rhs.push_back(nt->_actual_rhs[i]);
        v.push_back(nt->_actual_v[i]);
        area.push_back(nt->_actual_area[i]);

        if (tree.at(i).size()!=1) break; //end of branch of bifurcation
        i = tree.at(i).front(); //child below
    }
    a.shrink_to_fit();
    b.shrink_to_fit();
    d.shrink_to_fit();
    rhs.shrink_to_fit();
    v.shrink_to_fit();
    area.shrink_to_fit();

    //get all compartments information
    map<int,int> fromMechToIndex; //map from general mech ids to sequential ids
    vector<int> is_art, layout, is_ion;
    vector<char> pnttype;
    vector<Memb_func*> membfunc;
    vector<int> nodesIds, instanceCount, dataSize, pdataSize;
    vector<double> data;
    vector<Datum> pdata;

    int mechCount=0;
    for (map<int, list< pair<int, Memb_list*> > >::iterator mechs_it = mechanisms.begin(); mechs_it != mechanisms.end(); mechs_it++)
    {
        //per mechanism data
        int mechId = mechs_it->first();
        fromMechToIndex[mechId]=mechCount++;
        is_art.push_back(nrn_is_artificial_[mechId]);
        layout.push_back(nrn_mech_data_layout_[mechId]);
        pnttype.push_back(mechId);
        is_ion.push_back(nrn_is_ion(mechId));
        dataSize.push_back(nrn_prop_param_size_[mechId]);
        iDataSize.push_back(nrn_prop_dparam_size_[mechId]);
        membfunc.push_back(&memb_func[mechId]);

        //data related to the application of a mech to a node
        Memb_list* ml = mechs_it->second().second();
        int nodeId = fromNodeToIndex[mechs_it->second().first()];
        nodesIds.push_back(nodeId);
        instanceCount.push_back(ml->nodecount);

        data.insert (data.begin(),  ml->data,  ml->data+nrn_prop_param_size_[mechId]  );
        pdata.insert(pdata.begin(), ml->pdata, ml->pdata+nrn_prop_dparam_size_[mechId]);
    }
    assert(pnttype.size() == mechCount);
    is_art.shrink_to_fit();
    layout.shrink_to_fit();
    pnttype.shrink_to_fit();
    is_ion.shrink_to_fit();
    dataSize.shrink_to_fit();
    idataSize.shrink_to_fit();
    membfunc.shrink_to_fit();
    nodesIds.shrink_to_fit();
    instanceCount.shrink_to_fit();
    data.shrink_to_fit();
    pdata.shrink_to_fit();

    //bottom of branch = more than 1 children. Start recursive loop
    hpx_t * children = new hpx_t[nodes_it->second.length()];

    int c=0;
    for (list<int>::iterator it = tree.at(topNodeId).begin(); it != tre.at(topNodeId).end(); it++)
        children[c++] = nrx_setup_branch(nt, tree, *it);

    //end of branch, create it
    hpx_t branch_addr = hpx_gas_calloc_local(1, sizeof(Branch), NEUROX_HPX_MEM_ALIGNMENT);
    Branch branch;
    branch.children = children;
    branch.childrenCount = c;
    branch.n = n;
    branch.a = a.data();
    branch.b = b.data();
    branch.d = d.data();
    branch.area = area.data();
    branch.v = v.data();

    branch.is_art = is_art.data();
    branch.layout = layout.data();
    branch.is_ion = is_ion.data();
    branch.pnttype = pnttype.data();
    branch.membfunc = membfunc.data(); //todo not correctly serializable
    branch.nodesIds = nodesIds.data();
    branch.instanceCount = instanceCount.data();
    branch.dataSize = dataSize.data();
    branch.pdataSize = pdataSize.data();
    branch.data = data.data();
    branch.pdata = pdata.data();

    //TODO: send serialized version
    byte * bytes; int size;
    branch.serialize(bytes, size);
    hpx_call_sync(branch_addr, initializeBranch, NULL, 0, bytes, size);

    //return hpx_t of this branch
    delete [] children;
    return branch_addr;
}

void nrx_setup_neuron(NrnThread * nt, int gid, set<int> & neuronIds, hpx_t neuron_addr)
{
    //create tree of dependencies (key of top node is gid)
    map<int, list<int> > tree;
    for (set<int>::iterator it = neuronIds.begin(); it != neuronIds.end(); it++)
    {
        int i = *it;
        int p = nt->_v_parent_index[i];
        tree[p].push_back(i);
    }

    //create map of mechanisms (for a mech id, return the node id its memb_list)
    map<int, list< pair<int, Memb_list*> > > mechanisms;
    for (NrnThreadMembList* tml = nt->tml; tml; tml = tml->next)
    {
        int mechId = tml->index;
        Memb_list* ml = tml->ml; //mechanism's member list
        //if (!is_art) //TODO: is this valid?
        {
            for (int n=0; n< ml->nodecount; n++)
            {
                int nodeId = ml->nodeindices[n];
                mechanisms[mechId].push_back(std::make_pair(nodeId,ml));
            }
        }
    }

    //set-up neuron
    Neuron neuron;
    neuron.topBranch = nrx_setup_branch(nt, tree, mechanisms, gid);
    neuron.id = gid; //TODO should be global not local GID
    neuron.branches = children;

    //TODO serialization and GAS for neuron is done here
    hpx_call_sync(neuron_addr, initializeNeuron, done, &neuron, sizeof (Neuron));

    delete [] children;

    /*
    n=0;
    for (int k=0; k<nrn_nthread; k++)
    {
        NrnThread * nt = &nrn_threads[k];
        for (int i=nt->ncell; i<nt->end; i++)
        {
                int synOffset=0
                //point process for this mech
                if (pntype > 0)
                { //if is point proc
                   Point_process * pntprocs = nt.pntprocs[synOffset];
                   for (int s=0; s < ml->nodecount ; ++s)
                   {
                       Point_process *pntproc = &pntprocs[s];
                        //TODO: nt->_vdata is set as a pointer to this, used by mechs
                       pntproc->;
                   }
                   synOffset += ml->nodecount;
                }

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
        */

}

//all nodes have other data
void nrx_setup()
{
    unsigned firstNeuronGId=0; //gid of first neuron in this compute node

    //1. get total neurons count
    int neuronsCount=0;
    for (int k=0; k<nrn_nthread; k++)
    {
        FILE *dotfile = fopen(string("graph"+to_string(HPX_LOCALITY_ID)+"_"+to_string(k)+".dot").c_str(), "wt");
        fprintf(dotfile, "graph G%d\n{\n", k );

        //for all nodes in this NrnThread
        NrnThread * nt = &nrn_threads[k];
        neuronsCount += nt->ncell;
        for (int i=nt->ncell; i<nt->end; i++)
            fprintf(dotfile, "%d -- %d;\n", nt->_v_parent_index[i], i);
        fprintf(dotfile, "}\n");
        fclose(dotfile);
    }

    //TODO here they should reduce-all and get count of all nodes

    GlobalInfo info;
    info.neuronsCount = n; //TODO global count instead
    info.neuronsAddr  = hpx_gas_calloc_blocked(info.neuronsCount, sizeof(Neuron), NEUROX_HPX_MEM_ALIGNMENT);
    info.multiSplit   = HPX_LOCALITIES > info.neuronsCount;
    assert(info.neuronsAddr != HPX_NULL);

    //Broadcast input arguments
    printf("Broadcasting global info...\n");
    int e = hpx_bcast_rsync(initializeGlobalInfo, &info, sizeof (GlobalInfo));
    assert(e == HPX_SUCCESS);

    int neuronId=0;
    for (int k=0; k<nrn_nthread; k++) //for all set of neurons
    {
        NrnThread * nt = &nrn_threads[k];
        for (int n=0; n<nt->ncell; n++, neuronId++) //for all neurons in this NrnThread
        {
            //collect list of nodes from this neuron
            std::set<int> nodesFromThisNeuron;

            //insert soma Id
            nodesFromThisNeuron.insert(n);

            for (int i=nt->ncell; i<nt->end; i++)
                //if father belongs to this neuron, so do I
                if (nodesFromThisNeuron.find(nt->_v_parent_index[i]) != nodesFromThisNeuron.end())
                    nodesFromThisNeuron.insert(i);

            hpx_t neuron_addr = hpx_addr_add(globalInfo->neuronsAddr, sizeof(Neuron)*(neuronId+firstNeuronGid), sizeof(Neuron));
            nrx_setup_neuron(nt, n, nodesFromThisNeuron, neuron_addr);
        }
    }
}



void nrx_setup_register_hpx_actions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  initializeGlobalInfo, initializeGlobalInfo_handler, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  initializeNeuron, initianizeNeuron_handler, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  initializeBranch, initializeBranch_handler, HPX_POINTER, HPX_SIZE_T);
}


