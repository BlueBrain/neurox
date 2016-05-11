#include <stdio.h>
#include <string>
#include <set>
#include <map>
#include <list>
#include <tuple>

#include "coreneuron/nrnconf.h"
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrnoc/nrnoc_decl.h"
#include "coreneuron/nrniv/netcon.h"
#include "coreneuron/nrniv/netcvode.h"
#include "coreneuron/nrniv/nrniv_decl.h"
#include "coreneuron/nrniv/nrn_assert.h"

#include "neurox/neurox.h"

InputParams * inputParams;
Circuit * circuit;

using namespace std;

NrxSetup::NrxSetup()
{}

NrxSetup::~NrxSetup()
{}


hpx_t NrxSetup::createBranch(NrnThread * nt, map<int, list<int> > & tree, map<int, list < std::tuple< int, double*, int*> > > & mechanisms, int topNodeId )
{
    vector<double> a,b,d,rhs,v,area; //nodes information
    map<int,int> fromNodeToIndex; //map from general nodes ids to sequential ids

    int n=0, currentNode=topNodeId;
    while (1)
    {
        fromNodeToIndex[currentNode]=n++;

        //node level info
        a.push_back(nt->_actual_a[currentNode]);
        b.push_back(nt->_actual_b[currentNode]);
        d.push_back(nt->_actual_d[currentNode]);
        rhs.push_back(nt->_actual_rhs[currentNode]);
        v.push_back(nt->_actual_v[currentNode]);
        area.push_back(nt->_actual_area[currentNode]);

        if (tree.at(currentNode).size()!=1) break; //end of branch of bifurcation
        currentNode = tree.at(currentNode).front(); //child below
    }
    a.shrink_to_fit();
    b.shrink_to_fit();
    d.shrink_to_fit();
    rhs.shrink_to_fit();
    v.shrink_to_fit();
    area.shrink_to_fit();

    //get all compartments information
    vector<int> is_art, layout, is_ion;
    vector<char> pnttype;
    vector<Memb_func*> membfunc;
    vector<int> nodesIds, instancesCount, dataSize, pdataSize;
    vector<double> data;
    vector<Datum> pdata;

    for (map<int, list < std::tuple< int, double*, int*> > >::iterator mechs_it = mechanisms.begin(); mechs_it != mechanisms.end(); mechs_it++)
    {
        //per mechanism data
        int mechId = mechs_it->first;
        is_art.push_back(nrn_is_artificial_[mechId]);
        layout.push_back(nrn_mech_data_layout_[mechId]);
        pnttype.push_back(mechId);
        is_ion.push_back(nrn_is_ion(mechId));
        dataSize.push_back(nrn_prop_param_size_[mechId]);
        pdataSize.push_back(nrn_prop_dparam_size_[mechId]);
        membfunc.push_back(&memb_func[mechId]); //TODO is this enough to serialize it

        //data related to the list of instances of mechanism
        list< std::tuple< int, double*, int* > > & mechsInstanceList = mechs_it->second;
        instancesCount.push_back(mechsInstanceList.size());
        for ( list < tuple< int, double*, int* > >::iterator instance_it = mechsInstanceList.begin(); instance_it != mechsInstanceList.end(); instance_it++)
        {
            tuple < int, double*, int* > & instance = *instance_it;
            int nodeId = fromNodeToIndex[std::get<0>(instance)];
            nodesIds.push_back(nodeId);
            data.insert (data.begin(),  std::get<1>(instance), std::get<1>(instance)+nrn_prop_param_size_[mechId] );
            pdata.insert(pdata.begin(), std::get<2>(instance), std::get<2>(instance)+nrn_prop_dparam_size_[mechId]);
        }
    }
    is_art.shrink_to_fit();
    layout.shrink_to_fit();
    pnttype.shrink_to_fit();
    is_ion.shrink_to_fit();
    dataSize.shrink_to_fit();
    pdataSize.shrink_to_fit();
    membfunc.shrink_to_fit();
    nodesIds.shrink_to_fit();
    instancesCount.shrink_to_fit();
    data.shrink_to_fit();
    pdata.shrink_to_fit();

    //bottom of branch has more than 1 children. So we start recursive loop
    hpx_t * children = new hpx_t[tree.at(currentNode).size()];

    int childrenCount=0;
    for (list<int>::iterator children_it = tree.at(currentNode).begin(); children_it != tree.at(currentNode).end(); children_it++)
        children[childrenCount++] = createBranch(nt, tree, mechanisms, *children_it);

    //end of branch, create it
    hpx_t branch_addr = hpx_gas_calloc_local(1, sizeof(Branch), NEUROX_HPX_MEM_ALIGNMENT);
    Branch branch;
    branch.children = children;
    branch.childrenCount = childrenCount;
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
    branch.membfunc = membfunc.data(); //TODO not correctly serializable
    branch.nodesIds = nodesIds.data();
    branch.instanceCount = instancesCount.data();
    branch.dataSize = dataSize.data();
    branch.pdataSize = pdataSize.data();
    branch.data = data.data();
    branch.pdata = pdata.data();

    /* Serialize SYNAPSES
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

    //send serializable version and constructs Branch from it
    byte * bytes; int size;
    branch.serialize(bytes, size);
    hpx_call_sync(branch_addr, Branch::initialize, NULL, 0, bytes, size);

    //return hpx_t of this branch
    delete [] children;
    return branch_addr;
}

void NrxSetup::createNeuron(NrnThread * nt, int gid, set<int> & neuronIds, hpx_t neuron_addr)
{
    //create tree of dependencies (key of top node is gid)
    map<int, list<int> > tree;
    for (set<int>::iterator it = neuronIds.begin(); it != neuronIds.end(); it++)
    {
        int i = *it;
        int p = nt->_v_parent_index[i];
        tree[p].push_back(i);
    }

    //create map of mechanisms (for a mech id, return the list of its instances)
    map<int, list < std::tuple< int, double*, int*> > > mechanisms;
    for (NrnThreadMembList* tml = nt->tml; tml; tml = tml->next)
    {
        int mechId = tml->index;
        Memb_list* ml = tml->ml; //mechanism's member list
        int pdataOffset=0, dataOffset=0;
        //if (!is_art) //TODO: is this valid?
        {
            for (int n=0; n< ml->nodecount; n++)
            {
                int nodeId = ml->nodeindices[n];
                if (neuronIds.find(nodeId) != neuronIds.end()) // part of this neuron
                    mechanisms[mechId].push_back(std::make_tuple(
                        nodeId, &ml->data[dataOffset], &ml->pdata[dataOffset]
                    ));
                dataOffset  += nrn_prop_param_size_[mechId];
                pdataOffset += nrn_prop_dparam_size_[mechId];
            }
        }
    }

    //set-up neuron
    Neuron neuron;
    neuron.topBranch = createBranch(nt, tree, mechanisms, gid);
    neuron.id = gid; //TODO should be global not local GID

    //allocate Neuron in GAS
    hpx_call_sync(neuron_addr, Neuron::initialize, NULL, 0, &neuron, sizeof (Neuron));
}

//all nodes have other data
void NrxSetup::copyFromCoreneuronToHpx()
{
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

    Circuit circuit;
    circuit.neuronsCount = neuronsCount; //TODO global count instead
    circuit.neuronsAddr  = hpx_gas_calloc_blocked(circuit.neuronsCount, sizeof(Neuron), NEUROX_HPX_MEM_ALIGNMENT);
    circuit.multiSplit   = HPX_LOCALITIES > circuit.neuronsCount;
    assert(circuit.neuronsAddr != HPX_NULL);

    //Broadcast input arguments
    printf("Broadcasting Circuit info...\n");
    int e = hpx_bcast_rsync(Circuit::initialize, &circuit, sizeof (Circuit));
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

            hpx_t neuron_addr = hpx_addr_add(circuit.neuronsAddr, sizeof(Neuron)*neuronId, sizeof(Neuron));
            createNeuron(nt, n, nodesFromThisNeuron, neuron_addr);
        }
    }
}
