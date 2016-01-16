#include "coreneuron/hpx/data_converter.h"

#include "coreneuron/nrnconf.h"
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrnoc/nrnoc_decl.h"
#include "coreneuron/nrniv/netcon.h"
#include "coreneuron/nrniv/netcvode.h"
#include "coreneuron/nrniv/nrniv_decl.h"
#include "coreneuron/nrniv/nrn_assert.h"


void convert_from_coreneuron_to_hpx_datatypes()
{
    //map of parent node hpx address, per id (for branching)
    std::map<int, hpx_t> parentsAddrPerNode;
    std::map<int, std::deque<Memb_list*>> mechsPerNode;
    //std::map<int, std::queue<NrnThreadMembList*>> dependenciesPerMech;

    //for all NrnThread in this group
    for (int k=0; k<nrn_nthread; k++)
    {
        NrnThread * nt = &nrn_threads[k];
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

    }
}

