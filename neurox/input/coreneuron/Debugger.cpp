#include <stdio.h>
#include <string>
#include <set>
#include <map>
#include <list>
#include <tuple>
#include <algorithm>
#include <numeric>

#include "coreneuron/nrnconf.h"
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrnoc/nrnoc_decl.h"
#include "coreneuron/nrniv/netcon.h"
#include "coreneuron/nrniv/ivocvect.h"
#include "coreneuron/nrniv/netcvode.h"
#include "coreneuron/nrniv/nrniv_decl.h"
#include "coreneuron/nrniv/nrn_assert.h"
#include "coreneuron/nrniv/nrn_setup.h"
#include "coreneuron/utils/memory_utils.h"
#include "coreneuron/nrnoc/nrnoc_decl.h" //nrn_is_ion()
#include "coreneuron/nrniv/vrecitem.h"

#include "neurox/neurox.h"

#define EPSILON_EQ_FLOAT 0.000001

using namespace std;
using namespace neurox::Input;
using namespace neurox::Input::Coreneuron;

bool Debugger::isEqual(floble_t a, floble_t b, bool roughlyEqual)
{
    return roughlyEqual ? fabs(a-b) < EPSILON_EQ_FLOAT : a==b;
}

void Debugger::compareMechanismsFunctionPointers( std::list<NrnThreadMembList*> & uniqueMechs)
{
    printf("NDEBUG::comparing Mechanisms functions...\n");
    for (auto & tml : uniqueMechs)
    {
        Memb_func & mf_cn = memb_func[tml->index]; //coreneuron
        Memb_func & mf_nx = mechanisms[mechanismsMap[tml->index]]->membFunc; //neurox
        if (tml->index != CAP)
        { //we call MechFunctions::jacobCapacitance and currentCapacitance
          //so value is set. Coreneuron sets to null and calls nrn_cur_capacitance instead
          assert (mf_cn.jacob == mf_nx.jacob);
          assert (mf_cn.current == mf_nx.current);
        }
        //assert (mf_cn.alloc == mf_nx.alloc); //TODO never used?
        assert (mf_cn.destructor == mf_nx.destructor);
        assert (mf_cn.initialize == mf_nx.initialize);
        assert (mf_cn.state == mf_nx.state);
        assert (mf_cn.thread_cleanup_ == mf_nx.thread_cleanup_);
        assert (mf_cn.thread_mem_init_ == mf_nx.thread_mem_init_);
        assert (mf_cn.thread_table_check_ == mf_nx.thread_table_check_);
    }
}

void Debugger::fixed_step_minimal()
{
    nrn_fixed_step_minimal(); //from fadvance_core.c
    //nrn_fixed_step_group_minimal(); //from fadvance_core.c
    //both call fadvance_core.c::nrn_fixed_step_thread(NrnThread*)
}

void Debugger::compareBranch2(Branch * branch)
{
    assert(branch->soma); //only non-branched neurons
    int nrnThreadId = branch->nt->id;
    assert(sizeof(floble_t) == sizeof(double)); //only works with doubles!
    NrnThread & nt = nrn_threads[nrnThreadId];

    bool multiMex = branch->mechsGraph != NULL;

    assert(branch->nt->_t == nt._t);
    assert(branch->soma->threshold   == nt.presyns[0].threshold_);
    assert(branch->soma->thvar_index == nt.presyns[0].thvar_index_);
    assert(branch->soma->gid == nt.presyns[0].gid_);

    //vecplay
    assert(branch->nt->n_vecplay == nt.n_vecplay);
    for (int i=0; i<nt.n_vecplay;i++)
    {
        void * vpx_ptr = branch->nt->_vecplay[i];
        VecPlayContinuousX * vpx = reinterpret_cast<VecPlayContinuousX*>(vpx_ptr);

        PlayRecord * prc = (PlayRecord*) nt._vecplay[i];
        VecPlayContinuous * vpc = (VecPlayContinuous*) prc;

        assert(prc->ith_ == branch->nt->id); //single neuron per NrnThread
        assert(vpx->last_index_ == vpc->last_index_);
        assert(vpx->discon_index_ == vpc->discon_index_);
        assert(vpx->size_  == vpc->y_->size());
        assert(vpx->size_  == vpc->t_->size());
        assert(isEqual(*(vpx->pd_), *(vpc->pd_), multiMex));
        assert((vpx->discon_indices_==NULL) == (vpc->discon_indices_==NULL));
        for (size_t j=0; j<vpx->size_; j++)
        {
            assert(vpx->y_[j] == vpc->y_->data()[j]); //constants
            assert(vpx->t_[j] == vpc->t_->data()[j]); //constants
            if (vpx->discon_indices_ != NULL)
              {assert(vpx->discon_indices_[j] == vpc->discon_indices_->data()[j]);} //constants
        }
    }

    //make sure morphology is correct
    for (int i=0; i<6*branch->nt->end; i++)
    {   assert(isEqual(nt._data[i], branch->nt->_data[i], multiMex)); }

    //make sure mechs data is correct
    for (int i=0; i<nt._ndata; i++)
    {   assert(isEqual(nt._data[i], branch->nt->_data[i], multiMex)); }

    for (offset_t i=0; i<branch->nt->end; i++)
    {
        assert(nt._actual_a[i] == branch->nt->_actual_a[i]); //constants
        assert(nt._actual_b[i] == branch->nt->_actual_b[i]); //constants
        assert(nt._actual_area[i] == branch->nt->_actual_area[i]); //constants
        assert(isEqual(nt._actual_d[i], branch->nt->_actual_d[i], multiMex));
        assert(isEqual(nt._actual_v[i], branch->nt->_actual_v[i], multiMex));
        assert(isEqual(nt._actual_rhs[i], branch->nt->_actual_rhs[i], multiMex));
        if (branch->nt->_v_parent_index)
        {    assert(nt._v_parent_index[i] == branch->nt->_v_parent_index[i]); }
    }

    //make sure weights are correct
    assert(nt.n_weight == branch->nt->n_weight);
    //for (int n=0; n< nt.n_weight; n++)
    //{   assert(nt.weights[n] == branch->nt->weights[n]); } //order is changed!


    //make sure netcons are correct
    size_t netconsCount =0;
    for (auto nc_it : branch->netcons)
    {
        int srcGid = nc_it.first;
        netconsCount += nc_it.second.size();
    }
    assert(nt.n_netcon == netconsCount);

    //make sure number of output synapses per neuron is correct;
    for (int n = 0; n< nt.n_netcon; n++)
    {
        //TODO how to compare?
        //NetConX * nx = branch->net
        //assert(nt.netcons[n].active_ == branch->nt->netcons[n].active_);
        //assert(nt.netcons[n].delay_ == branch->nt->netcons[n].delay_);
        //assert(nt.netcons[n].u.weight_index_ == branch->nt->netcons[n].delay_);
        //assert(nt.netcons[n].target_->_i_instance == branch->nt->netcons[n].delay_);
        //assert(nt.netcons[n].target_->_type == branch->nt->netcons[n].delay_);
    }

    int vdataOffset=0;
    for (NrnThreadMembList* tml = nt.tml; tml!=NULL; tml = tml->next) //For every mechanism
    {
        int type = tml->index;
        int m = mechanismsMap[type];
        Memb_list * ml = tml->ml; //Mechanisms application to each compartment
        Memb_list & instances = branch->mechsInstances[m];
        assert(ml->nodecount == instances.nodecount);
        //assert(ml->_nodecount_padded == instance.instancesCount);
        short dataSize  = mechanisms[m]->dataSize;
        short pdataSize = mechanisms[m]->pdataSize;
        for (int n=0; n<ml->nodecount; n++) //for every mech instance
        {
            assert(ml->nodeindices[n]==instances.nodeindices[n]);
            for (int i=0; i<dataSize; i++)
            {   assert(isEqual(ml->data[n*dataSize+i], instances.data[n*dataSize+i], multiMex)); }

            for (int i=0; i<pdataSize; i++)
            {
                int ptype = memb_func[type].dparam_semantics[i];
                assert(isEqual(ml->pdata[n*pdataSize+i], instances.pdata[n*pdataSize+i], multiMex));
            }

            /* We comment this because it runs for NULL presyn
            if (mechanisms[m]->pntMap)
            {
                //compare point_processes (index 1)
                Point_process * pp = (Point_process *) nt._vdata[vdataOffset+1];
                Point_process * pp2 = (Point_process *) branch->vdata[vdataOffset+1];
                assert(pp->_type == pp2->_type );
                assert(pp->_i_instance == pp2->_i_instance );
                assert(pp->_tid == pp2->_tid );
                assert(pp->_presyn == pp2->_presyn );
                vdataOffset+= mechanisms[m]->vdataSize;
            }
            */
        }
    }
}

hpx_action_t Debugger::compareBranch = 0;
int Debugger::compareBranch_handler()
{
    neurox_hpx_pin(Branch);
    compareBranch2(local);
    neurox_hpx_unpin;
}

void Debugger::coreNeuronFinitialize()
{
    nrn_finitialize(inputParams->voltage != 1000., inputParams->voltage);
}

void Debugger::registerHpxActions()
{
    neurox_hpx_register_action(0, Debugger::compareBranch);
}
