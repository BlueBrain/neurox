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

using namespace std;
using namespace neurox::input;

bool Debugger::IsEqual(floble_t a, floble_t b, bool roughlyEqual)
{
    const double epsilon = inputParams->multiMex ? 5e-4 : 1e-8;
    //NOTE: current functions for mechs SK_E2 (144) andNap_Et2 (123) depends on
    //V that change opening vars that change V and so on... after few steps this
    //causes high difference in results;
    return roughlyEqual ? fabs(a-b) < epsilon : a==b;
}

void Debugger::CompareMechanismsFunctionPointers()
{
#if !defined(NDEBUG)
    printf("NDEBUG::comparing Mechanisms functions...\n");
    for (int m=0; m<neurox::mechanismsCount; m++)
    {
        int index = mechanisms[m]->type;
        Memb_func & mf_cn = memb_func[index]; //coreneuron
        Memb_func & mf_nx = mechanisms[m]->membFunc; //neurox

        //constructur, destructors are different
        //assert (mf_cn.alloc == mf_nx.alloc);
        //assert (mf_cn.destructor == mf_nx.destructor);

        //file reader is not used (we read directly from Coreneuron)
        //assert (mechanisms[m]->nrn_bbcore_read == nrn_bbcore_read_[index]);

        //cap calls MechFunctions::jacobCapacitance and currentCapacitance instead
        if (index != CAP)
        {
          assert (mf_cn.current == mf_nx.current);
          assert (mf_cn.jacob == mf_nx.jacob);
        }

        assert (mf_cn.state == mf_nx.state);
        assert (mf_cn.initialize == mf_nx.initialize);

        //If fails below, probably forgot to disable threading table when generating circuit
        assert (mf_cn.setdata_ == mf_nx.setdata_);
        assert (mf_cn.thread_cleanup_ == mf_nx.thread_cleanup_);
        assert (mf_cn.thread_mem_init_ == mf_nx.thread_mem_init_);
        assert (mf_cn.thread_size_ == mf_nx.thread_size_);
        assert (mf_cn.thread_table_check_ == mf_nx.thread_table_check_);

        for (int i=0; i< BEFORE_AFTER_SIZE; i++)
            if (nrn_threads[0].tbl[i])
                assert (mechanisms[m]->BAfunctions[i] ==  nrn_threads[0].tbl[i]->bam->f);
    }
#endif
}

void Debugger::StepAfterStepFinitialize(Branch *b, NrnThread *nth)
{
    b->CallModFunction(Mechanism::ModFunction::threadTableCheck);
    b->InitVecPlayContinous();

    //coreneuron
    t = 0.;
    dt2thread(-1.);
    nrn_thread_table_check();
    clear_event_queue();
    nrn_spike_exchange_init();
    nrn_play_init();

    /****************/ CompareBranch2(b); /*****************/

    b->DeliverEvents(b->nt->_t);
    for (int n=0; n<b->nt->end; n++)
        b->nt->_actual_v[n]=inputParams->voltage;

    //coreneuron
    nrn_deliver_events(nth); /* The play events at t=0 */

    for (int i = 0; i < nth->end; ++i)
        nth->_actual_v[i] = inputParams->voltage;

    /****************/ CompareBranch2(b); /*****************/

    b->CallModFunction(Mechanism::ModFunction::before_initialize);

    nrn_ba(nth, BEFORE_INITIAL);

    /****************/ CompareBranch2(b); /*****************/

    b->CallModFunction(Mechanism::ModFunction::initialize);

    NrnThreadMembList* tml;
    for (tml = nth->tml; tml; tml = tml->next) {
        mod_f_t s = memb_func[tml->index].initialize;
        if (s) {
            (*s)(nth, tml->ml, tml->index);
            }
        }
    init_net_events();

    /****************/ CompareBranch2(b); /*****************/

    b->CallModFunction(Mechanism::ModFunction::after_initialize);
    b->DeliverEvents(b->nt->_t);

    //coreneuron

    nrn_ba(nth, AFTER_INITIAL);
    nrn_deliver_events(nth); /* The INITIAL sent events at t=0 */

    /****************/ CompareBranch2(b); /*****************/

    b->SetupTreeMatrix();

    setup_tree_matrix_minimal(nth);

    /****************/ CompareBranch2(b); /*****************/

    b->CallModFunction(Mechanism::ModFunction::before_step);
    b->DeliverEvents(b->nt->_t);

    nrn_ba(nth, BEFORE_STEP);
    nrn_deliver_events(nth);

    /****************/ CompareBranch2(b); /*****************/
}

void Debugger::StepAfterStepBackwardEuler(Branch *b, NrnThread * nth, int secondorder)
{
    double dt = b->nt->_dt;
    if (b->soma && inputParams->algorithm==neurox::Algorithm::BackwardEulerWithTimeDependencyLCO)
    {
        b->soma->timeDependencies->SendSteppingNotification(b->nt->_t, dt, b->soma->gid, b->soma->synapses);
        b->soma->timeDependencies->WaitForTimeDependencyNeurons(b->nt->_t, dt, b->soma->gid);
    }
    if (b->soma)
    {
      //Soma waits for AIS to have threshold V value updated
      floble_t thresholdV;
      solver::HinesSolver::SynchronizeThresholdV(b, &thresholdV);
      if (b->soma->CheckAPthresholdAndTransmissionFlag(thresholdV))
          b->soma->SendSpikes(b->nt->_t);
          //TODO sendSpikes LCO must be waited
    }
    else if (b->thvar_ptr)
        //Axon Initial Segment send threshold  V to parent
        solver::HinesSolver::SynchronizeThresholdV(b);

    b->nt->_t += .5*dt;
    b->DeliverEvents(b->nt->_t);

    //coreneuron
    assert(b->nt->cj == nth->cj);
    dt2thread(dt);;
    deliver_net_events(nth);
    nth->_t += .5 * nth->_dt;

    /****************/ CompareBranch2(b); /*****************/

    b->FixedPlayContinuous();
    b->SetupTreeMatrix();
    b->SolveTreeMatrix();

    //coreneuron
    fixed_play_continuous(nth);
    setup_tree_matrix_minimal(nth);
    nrn_solve_minimal(nth);

    /****************/ CompareBranch2(b); /*****************/

    second_order_cur(b->nt, inputParams->secondorder );
    floble_t secondOrderMultiplier = inputParams->secondorder ? 2.0 : 1.0;
    for (int i=0; i<b->nt->end; i++)
        b->nt->_actual_v[i] += secondOrderMultiplier * b->nt->_actual_rhs[i];
    b->CallModFunction(Mechanism::ModFunction::currentCapacitance);

    second_order_cur(nth, secondorder);
    update(nth);

    /****************/ CompareBranch2(b); /*****************/

    b->nt->_t += .5*dt;
    b->FixedPlayContinuous();

    nth->_t += .5 * nth->_dt;//nrn_fixed_step_lastpart(nth);
    fixed_play_continuous(nth);

    /****************/ CompareBranch2(b); /*****************/

    b->CallModFunction(Mechanism::ModFunction::state);

    nonvint(nth);

    /****************/ CompareBranch2(b); /*****************/

    b->CallModFunction(Mechanism::ModFunction::after_solve);
    b->CallModFunction(Mechanism::ModFunction::before_step);
    b->DeliverEvents(b->nt->_t);

    nrn_ba(nth, AFTER_SOLVE);
    nrn_ba(nth, BEFORE_STEP);
    nrn_deliver_events(nth);

    /****************/ CompareBranch2(b); /*****************/
}

//fadvance_core.c::fixed_step_minimal(...)
void Debugger::FixedStepMinimal2(NrnThread * nth, int secondorder)
{
    //fadvance_core.c::dt2thread
    if (secondorder)
        nth->cj = 2.0 / nth->_dt;
    else
        nth->cj = 1.0 / nth->_dt;

    //fadvance_core.c::nrn_fixed_step_thread(...)
    deliver_net_events(nth);
    nth->_t += .5 * nth->_dt;
    if (nth->ncell)
    {
        fixed_play_continuous(nth);
        setup_tree_matrix_minimal(nth);
        nrn_solve_minimal(nth);
        second_order_cur(nth, secondorder);
        update(nth);
    }
    if (!nrn_have_gaps) {
        nrn_fixed_step_lastpart(nth);
    }
}

void Debugger::CompareAllBranches()
{
#if !defined(NDEBUG)
    if (inputParams->branchingDepth>0) return;
    DebugMessage("neurox::Input::CoreNeuron::Debugger::compareBranch...\n");
    neurox_hpx_call_neurons(input::Debugger::CompareBranch);
#endif
}

void Debugger::CompareBranch2(Branch * branch)
{
    assert(branch->soma); //only non-branched neurons
    int nrnThreadId = branch->nt->id;
    assert(sizeof(floble_t) == sizeof(double)); //only works with doubles!
    NrnThread & nt = nrn_threads[nrnThreadId];

    bool multiMex = branch->mechsGraph != NULL;

    assert(branch->nt->_t == nt._t);
    assert(branch->nt->_ndata == nt._ndata);
    assert(secondorder == inputParams->secondorder);
    assert(branch->soma->threshold   == nt.presyns[0].threshold_);
    assert(branch->soma->gid == nt.presyns[0].gid_);
    assert(IsEqual(*(branch->thvar_ptr), nt._actual_v[nt.presyns[0].thvar_index_], multiMex));

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
        assert(IsEqual(*(vpx->pd_), *(vpc->pd_), multiMex));
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
    assert(branch->nt->end == nt.end);
    for (int i=0; i<6; i++) //RHS, D, A, B, V and area
        for (int j=0; j<branch->nt->end; j++) //for all compartments
        {
            int offset = tools::Vectorizer::SizeOf(nt.end)*i+j;
            assert(IsEqual(nt._data[offset], branch->nt->_data[offset], multiMex));
        }

    for (offset_t i=0; i<branch->nt->end; i++)
    {
        assert(nt._actual_a[i] == branch->nt->_actual_a[i]); //constants
        assert(nt._actual_b[i] == branch->nt->_actual_b[i]); //constants
        assert(nt._actual_area[i] == branch->nt->_actual_area[i]); //constants
        assert(IsEqual(nt._actual_d[i], branch->nt->_actual_d[i], multiMex));
        assert(IsEqual(nt._actual_v[i], branch->nt->_actual_v[i], multiMex));
        assert(IsEqual(nt._actual_rhs[i], branch->nt->_actual_rhs[i], multiMex));
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
            {
#if LAYOUT==1
                int offset = mechanisms[m]->dataSize*n+i;
#else
                int offset = tools::Vectorizer::SizeOf(ml->nodecount)*i+n ;
#endif
                assert(IsEqual(ml->data[offset], instances.data[offset], multiMex));
            }

            for (int i=0; i<pdataSize; i++)
            {
#if LAYOUT==1
                int offset = pdataSize*n+i;
#else
                int offset = tools::Vectorizer::SizeOf(ml->nodecount)*i+n ;
#endif
                int ptype = memb_func[type].dparam_semantics[i];
                bool isPointer = ptype==-1 || (ptype>0 && ptype<1000);
                if (isPointer)
                {   assert(IsEqual(nt._data[ml->pdata[offset]], branch->nt->_data[instances.pdata[offset]], multiMex)); }
                assert(ml->pdata[offset] == instances.pdata[offset]);
            }

            /* We comment this because it runs for NULL presyn
            if (mechanisms[m]->pntMap)
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

hpx_action_t Debugger::FixedStepMinimal = 0;
int Debugger::FixedStepMinimal_handler(const int *steps_ptr, const size_t)
{
    neurox_hpx_pin(uint64_t);
    for (int n=0; n < nrn_nthread; n++)
      for (int i=0; i< *steps_ptr; i++)
        Debugger::FixedStepMinimal2(&nrn_threads[n], inputParams->secondorder);
    neurox_hpx_unpin;
}

hpx_action_t Debugger::Finitialize = 0;
int Debugger::Finitialize_handler()
{
    neurox_hpx_pin(uint64_t);
    nrn_finitialize(inputParams->voltage != 1000., inputParams->voltage);
    neurox_hpx_unpin;
}

hpx_action_t Debugger::ThreadTableCheck = 0;
int Debugger::ThreadTableCheck_handler()
{
    //beginning of fadvance_core.c::nrn_fixed_step_group_minimal
    neurox_hpx_pin(uint64_t);
    dt2thread(dt);  //does nothing
    nrn_thread_table_check();
    neurox_hpx_unpin;
}

hpx_action_t Debugger::NrnSpikeExchange = 0;
int Debugger::NrnSpikeExchange_handler()
{
    neurox_hpx_pin(uint64_t);
    nrn_spike_exchange(nrn_threads);
    neurox_hpx_unpin;
}

hpx_action_t Debugger::CompareBranch = 0;
int Debugger::CompareBranch_handler()
{
    neurox_hpx_pin(Branch);
    if (inputParams->branchingDepth==0)
        CompareBranch2(local); //not implemented for branch-parallelism
    neurox_hpx_unpin;
}

void Debugger::RegisterHpxActions()
{
    neurox_hpx_register_action(neurox_zero_var_action,   Debugger::CompareBranch);
    neurox_hpx_register_action(neurox_zero_var_action,   Debugger::Finitialize);
    neurox_hpx_register_action(neurox_zero_var_action,   Debugger::NrnSpikeExchange);
    neurox_hpx_register_action(neurox_zero_var_action,   Debugger::ThreadTableCheck);
    neurox_hpx_register_action(neurox_single_var_action, Debugger::FixedStepMinimal);
}

