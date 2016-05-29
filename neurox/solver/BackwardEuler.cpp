#include "neurox/neurox.h"
#include <numeric>
#include <algorithm>

BackwardEuler::BackwardEuler(InputParams * inputParams)
    :t(0)
{
    dt=inputParams->dt;
    tstop->inputParams->tstop;
}

BackwardEuler::~BackwardEuler() {}

int BackwardEuler::solve()
{
    for (double t=0; t<tstop; t+=dt)
        hpx_par_for_sync( [&] (int i, void*)
            { hpx_call_sync(brain->getNeuronAddr(i), BackwardEuler::step, NULL, 0);},
            0, brain->neuronsCount,  dt, sizeof(double));
}

hpx_action_t BackwardEuler::step = 0;
int BackwardEuler::step_handler(const double dt)
{
    neurox_hpx_pin(Neuron); //We are in a Neuron

    //multicore.c::nrn_thread_table_check()
    local->callMechsFunction_handler(Mechanism::Function::tableCheck);

    //multicore.c::nrn_fixed_step_thread()

    //cvodestb::deliver_net_events(nth);
    //1. send outgoing spikes netcvode.cpp::NetCvode::check_thresh()
    hpx_addr_t synapsesLco = HPX_NULL;
    if (true) //TODO check APT:
    //if (local->topBranch->v[0] local->thresholdAP)
    {
        //netcvode.cpp::PreSyn::send()
        synapsesLco = hpx_lco_and_new(local->synapsesCount);
        hpx_par_for( [&] (int i, void*)
        { hpx_call(local->synapses[i], BackwardEuler::queueSpike,
                   synapsesLco, local->id, sizeof(int));},
        0, local->synapsesCount, NULL);
    }
    hpx_lco_wait(synapsesLco);
    hpx_lco_delete(synapsesLco, HPX_NULL);
    //TODO lco should be until all synapses are received, not sent!

    //2. netcvode.cpp::NetCvode::deliver_net_events()
    hpx_call_sync(local->topBranch, BackwardEuler::deliverSpikes, NULL, 0);

    local->t += .5*dt;

    //fixed_play_continuous; for PlayRecord

    //netpar.c::nrn_spike_exchange()

    //TODO: nrn_thread_table_check() missing;

    //process queued synapses
    //(*pnt_receive_t)(Point_process*, double*, double)

    //run step

    //send synapses


    neurox_hpx_unpin;
}


hpx_action_t BackwardEuler::queueSpike = 0;
int BackwardEuler::queueSpike_handler(const Synapse * syn, size_t)
{
    neurox_hpx_pin(Branch);
    //netcvode::PreSyn::send()
    //TODO:Coreneuron uses a priority queue, why? (tqueue.c)
    hpx_lco_sema_p(local->mutex);
    local->queuedSynapses.push(*syn);
    hpx_lco_sema_v_sync(local->mutex);
    neurox_hpx_unpin;
}

hpx_action_t BackwardEuler::deliverSpikes = 0;
int BackwardEuler::deliverSpikes_handler()
{
    neurox_hpx_pin(Branch);
    //netcvode.cpp::NetCvode::deliver_net_events()
    //            ->NetCvode::deliver_events()
    //            ->NetCvode::deliver_event()
    //            ->NetCon::deliver()

    //Launch in all sub branches
    hpx_addr_t lco = hpx_lco_and_new(local->childrenCount);
    for (int c=0; c<local->childrenCount; c++)
        hpx_call(local->children[c], BackwardEuler::deliverSpikes,lco);

    //Deliver all spikes
    while (local->queuedSynapses.size())
    {
        Synapse & syn = local->queuedSynapses.front();
        (*brain->mechsTypes[syn.mechType].pnt_receive_t)((void *) syn.mechInstance, (void *) syn.weight, syn.delay);
        local->queuedSynapses.pop();
    }
    hpx_lco_wait(lco);
    hpx_lco_delete(lco, HPX_NULL);
    neurox_hpx_unpin;
}

static void BackwardEuler::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE, solve, solve_handler);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED, step, step_handler, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED, queueSpike, queueSpike_handler, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED, deliverSpikes, deliverSpikes_handler);
}

