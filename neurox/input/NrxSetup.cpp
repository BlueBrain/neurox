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

#include "neurox/neurox.h"

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

void NrxSetup::createNeuron(NrnThread & nt, int nrnThreadId, int neuronId)
{
    //set-up neuron
    Neuron neuron;
    //neuron.topBranch = createBranch(nt, tree, mechanisms, gid);
    neuron.id = gid; //TODO should be global not local GID

    //allocate Neuron in GAS
    hpx_call_sync(neuron_addr, Neuron::initialize, NULL, 0, &neuron, sizeof (Neuron));
}

void NrxSetup::copyFromCoreneuronToHpx()
{
    Brain brain;

    //1. get total neurons count; plot neurons tree to dot file
    //std::for_each(nrn_threads, nrn_threads+nrn_nthread, [&](NrnThread & nt){ circuit.neuronsCount+= nt.ncell; });
    brain.neuronsCount = std::accumulate(nrn_threads, nrn_threads+nrn_nthread, 0, [](int n, NrnThread & nt){return n+nt.ncell;});
    brain.neuronsAddr = hpx_gas_calloc_blocked(brain.neuronsCount, sizeof(Neuron), NEUROX_HPX_MEM_ALIGNMENT);
    brain.multiSplit = HPX_LOCALITIES > brain.neuronsCount;
    assert(brain.neuronsAddr != HPX_NULL);

    //TODO: Debug: plot morpholog as dot file
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

    //Broadcast input arguments
    printf("Broadcasting Circuit info...\n");
    int e = hpx_bcast_rsync(Brain::initialize, &brain, sizeof (Brain));
    assert(e == HPX_SUCCESS);

    //create the tree structure for all neurons
    vector<shared_ptr<Compartment>> compartments;

    long long acc=0; //accumulator of compartments ids
    for_each(nrn_threads, nrn_threads+nrn_nthread, [&](NrnThread & nt, int & id) {

        //reconstructs tree with solver values
        //for (int n=nt.end-1; n>=nt.ncell; n--)
        for (int n=nt.end-1; n>=0; n--)
        {
            Compartment & compartment = compartments[n]+acc;
            compartment.setSolverValues(nt._actual_a[n], nt._actual_d[n], nt._actual_v[n], nt._actual_d[n], nt._actual_rhs[n], nt._actual_area[n]);
            //if it is not top node, i.e. has a parent
            if ( n>= nt.ncell )
            {
                Compartment & parentCompartment = compartment[nt._v_parent_index[n]+acc];
                parentCompartment.addChild(&compartment);
            }
        }

        //reconstructs mechanisms for compartments
        int dataOffset=0, pdataOffset=0;
        for (NrnThreadMembList* tml = nt.tml; tml!=nullptr; tml = tml->next)
        {
            int mechId = tml->index;
            if (nrn_is_artificial_[mechId])
            {
                //TODO
            }
            if (pnt_map[mechId])
            {
                //TODO
            }

            Memb_list *& ml = tml->ml; //list of compartments this mechanism belongs to
            int & dataSize  = nrn_prop_param_size_[mechId];
            int & pdataSize = nrn_prop_dparam_size_[mechId];
            for (int n=0; n<ml->nodecount; n++)
            {
                Compartment & compartment = compartments[ml->nodeindices[n]+acc];
                compartment.mechanismsIds.insert(tml->index);
                compartment.data.insert (compartment.data.end(),  &ml->data[dataOffset],  &ml->data[dataOffset+dataSize]);
                compartment.pdata.insert(compartment.pdata.end(), &ml->pdata[pdataOffse], &ml->pdata[pdataOffset+pdataSize]);
                dataOffset  += dataSize;
                pdataOffset += pdataSize;
            }
        }

        //increments accumulator so that next NrnThread ids don't overlapt current one
        //acc += nt.end - nt.ncell;
        acc += nt.ncell;
    });

    //Builds mechanisms dependencies graph: 1: count dependencies
    MechanismsDependencies mechs;
    NrnThread & nt = nrn_threads[0];
    for (NrnThreadMembList* tml = nt.tml; tml!=nullptr; tml = tml->next)
         mechs.count[tml->index]++;

    for (int i=0; i<= mechs.count.size(); i++) //for all mechanisms
    {
        mechs.offsets[i]= i==0 ? 0 : mechs.offsets[i-1] + mechs.count[i-1];
        offset += mechs.offsets[i];
    }

    createNeuron(nt, n, id); //TODO create and distribute neuron as we read
}
