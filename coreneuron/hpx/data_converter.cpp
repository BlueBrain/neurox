#include "coreneuron/hpx/data_converter.h"

#include "coreneuron/nrnconf.h"
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrnoc/nrnoc_decl.h"
#include "coreneuron/nrniv/netcon.h"
#include "coreneuron/nrniv/netcvode.h"
#include "coreneuron/nrniv/nrniv_decl.h"
#include "coreneuron/nrniv/nrn_assert.h"


#define PP2NT(pp) (nrn_threads + (pp)->_tid)

void convert_from_coreneuron_to_hpx_datatypes(double tt, NetCvode* ns, NrnThread* nt)
{
    //map of parent node hpx address, per id (for branching)
    std::map<int, hpx_t> parentsAddrPerNode;
    std::map<int, std::queue<Memb_list*>> mechsPerNode;
    //std::map<int, std::queue<NrnThreadMembList*>> dependenciesPerMech;

    //for all NrnThread in this group
    for (int t=0; t<nrn_nthread; t++)
    {
        NrnThread * nt = &nrn_threads[t];
        double * nt_data = nt->_data; //nt->_ndata;
        //nidata, etc etc

        //for all mechanisms on this NrnThread
        for (NrnThreadMembList* tml = nt.tml; tml; tml = tml->next)
        {
            int type = tml->index;

            int is_art = nrn_is_artificial_[type];
            int layout = nrn_mech_data_layout_[mtype];
            char pntype = pnt_map[type];
            int is_ion = nrn_is_ion(type); //mem_func[type]...

            //Member functions (see declaration of initialize, etc)
            Memb_func* memfunc = &memb_func[index];

            /* TODO: dependencies not used on the code */
            //for (int d=0; d<tml->ndependencies; d++)
            //    dependenciesPerMech[tml].push_back( tml->dependencies[d]);

            Memb_list* ml = tml->ml;
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

        //for all nodes in this NrnThread
        for (int i=nt->ncell; i<nt->end; i++)
        {
            //TODO: what's the p of the very first node
            int    p    = nt->_v_parent_index[i];
            double a    = nt->_actual_a[i];
            double b    = nt->_actual_b[i];
            double d    = nt->_actual_d[i];
            double rhs  = nt->_actual_rhs[i];
            double v    = nt->_actual_v[i];
            double area = nt->_actual_area[i];

            //TODO: reverse engineer morphologY
        }

        //traverse preSyns (local spike exchange / event delivery)
        PreSyn * nt->presyns;
        for (int s=0; s<nt->n_presyn; n++)
        {
            PreSyn * preSyn = nt->presyns[s];
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
        gid2

        nt->tbl; //table for mechanism updates
        nt->_ecell_memb_list; //??

    }


    PreSyn * a;
    for (int i = a->nc_cnt_-1; i >= 0; --i) {
        NetCon* d = netcon_in_presyn_order_[a->nc_index_ + i];
        if (d->active_ && d->target_) {
            NrnThread* n = PP2NT(d->target_);
            if (nt == n) {
                ns->bin_event(tt + d->delay_, d, n);
            }else{
                ns->p[n->id].interthread_send(tt + d->delay_, d, n);
            }
        }
        }
}

