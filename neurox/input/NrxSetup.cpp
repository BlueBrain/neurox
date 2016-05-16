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

void NrxSetup::copyFromCoreneuronToHpx()
{
    int neuronsCount =  std::accumulate(nrn_threads, nrn_threads+nrn_nthread, 0, [](int n, NrnThread & nt){return n+nt.ncell;});
    createBrain(neuronsCount);

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
    int acc=0;
    for_each(nrn_threads, nrn_threads+nrn_nthread, [&](NrnThread & nt, int & acc) {

        vector<Compartment> compartments;
        vector<Mechanism> mechanisms;

        //reconstructs tree with solver values
        for (int n=nt.end-1; n>=0; n--)
        {
            Compartment & compartment = compartments[n];
            compartment.setSolverValues(nt._actual_a[n], nt._actual_b[n], nt._actual_d[n], nt._actual_v[n], nt._actual_rhs[n], nt._actual_area[n]);

            if ( n>=nt.ncell) //if it is not top node, i.e. has a parent
            {
                Compartment & parentCompartment = compartment[nt._v_parent_index[n]];
                parentCompartment.addChild(&compartment);
            }
        }

        //reconstructs mechanisms for compartments
        int dataOffset=0, pdataOffset=0;
        for (NrnThreadMembList* tml = nt.tml; tml!=nullptr; tml = tml->next)
        {
            int mechId = tml->index;
            //Data unique to each mechanism type
            mechanisms[mechId](nrn_prop_param_size_[mechId],
                    nrn_prop_dparam_size_[mechId],tml->ndependencies,
                    pnt_map[mechId], nrn_is_artificial_[mechId], tml->dependencies);

            //Mechanisms' application to each compartment
            Memb_list *& ml = tml->ml; //list of compartments this mechanism belongs to
            int & dataSize  = mechanisms[mechId].dataSize;
            int & pdataSize = mechanisms[mechId].pdataSize;
            for (int n=0; n<ml->nodecount; n++)
            {
                Compartment & compartment = compartments[ml->nodeindices[n]];
                compartment.mechanismsIds.insert(mechId);
                compartment.data.insert (compartment.data.end(),  &ml->data[dataOffset],  &ml->data[dataOffset+dataSize]);
                compartment.pdata.insert(compartment.pdata.end(), &ml->pdata[pdataOffse], &ml->pdata[pdataOffset+pdataSize]);
                dataOffset  += dataSize;
                pdataOffset += pdataSize;
            }
        }

        //We will now convert from compartment level to branch level
        for (int n=0; n<nt.end; n++)
            createNeuron(acc+n, compartments[n], mechanisms);
        acc +=n;
    });
}

void NrxSetup::createBrain(int neuronsCount)
{
    //create Brain HPX structure
    Brain brain;
    brain.neuronsCount = neuronsCount;
    brain.neuronsAddr = hpx_gas_calloc_blocked(brain.neuronsCount, sizeof(Neuron), NEUROX_HPX_MEM_ALIGNMENT);
    brain.multiSplit = HPX_LOCALITIES > brain.neuronsCount;
    assert(brain.neuronsAddr != HPX_NULL);

    //Broadcast input arguments
    printf("Broadcasting Circuit info...\n");
    int e = hpx_bcast_rsync(Brain::initialize, &brain, sizeof (Brain));
    assert(e == HPX_SUCCESS);
}

void NrxSetup::createNeuron(int gid, vector<Compartment> & compartments, vector<Mechanism> & mechanisms)
{
    Neuron neuron;
    neuron.topBranch = createBranch(&compartments[n], mechanisms);
    neuron.id = gid;

    //allocate Neuron in GAS
    hpx_call_sync(brain[gid], Neuron::initialize, NULL, 0, &neuron, sizeof (Neuron));
}

hpx_t NrxSetup::createBranch( Compartment * topCompartment, vector<Mechanism> & mechanisms)
{
    Compartment *comp = topCompartment;
    Branch b;
    b.n =0; //number of compartments
    for (Compartment *comp = topCompartment;
         comp.right == nullptr && comp.left != nullptr;
         comp = comp->left, b.n++)
    {
        //TODO copy vars
    };

    b.childrenCount = (comp->right==nullptr ? 0 : 1) + (comp->left==nullptr ? 0 : 1);
    if (b.childrenCount > 0)
    {
        b.children = new hpx_t[2];
        b.children[0] = createBranch (comp->left);
        b.children[1] = createBranch (comp->right);
    }
    else b.children = nullptr;

    //Allocate HPX Branch
    hpx_t branchAddr = hpx_gas_calloc_local(1, sizeof(Branch), NEUROX_HPX_MEM_ALIGNMENT);
    hpx_call_sync(branchAddr, Branch::initialize, NULL, 0, bytes, size);
    return branchAddr;

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
}
