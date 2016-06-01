#include "neurox/Neurox.h"
#include <numeric>
#include <algorithm>

using namespace Neurox;
using namespace Neurox::Solver;

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

    //1. multicore.c::nrn_thread_table_check()
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::Function::tableCheck);

    //2. multicore.c::nrn_fixed_step_thread()

    //2.1 cvodestb::deliver_net_events(nth);
    //(send outgoing spikes netcvode.cpp::NetCvode::check_thresh() )
    //TODO check APT: if (local->topBranch->v[0] local->thresholdAP)
    bool reachedThresold = true;
    hpx_addr_t spikesLco = reachedThresold ? Neuron::fireActionPotential(local) : HPX_NULL;

    //3. netcvode.cpp::NetCon::deliver()
    //calls NET_RECEIVE in mod files to receive synapses
    //TODO: are these the same ones sent, or from the previous interval?
    hpx_call_sync(local->topBranch, BackwardEuler::deliverSpikes, NULL, 0);

    local->t += .5*dt;

    //TODO: fixed_play_continuous; for PlayRecord (whole logic missing)

    Neuron::setupTreeMatrixMinimal(local);

    //netpar.c::nrn_spike_exchange()

    //TODO: nrn_thread_table_check() missing;

    //process queued synapses
    //(*pnt_receive_t)(Point_process*, double*, double)

    //run step

    //send synapses


    //wait for all post-synaptic neurons to receive (not process) synapses
    if (spikesLco != HPX_NULL)
    {
        hpx_lco_wait(spikesLco);
        hpx_lco_delete(spikesLco, HPX_NULL);
    }
    neurox_hpx_unpin;
}





static void BackwardEuler::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE, solve, solve_handler);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED, step, step_handler, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED, queueSpike, queueSpike_handler, HPX_POINTER, HPX_SIZE_T);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED, deliverSpikes, deliverSpikes_handler);
}

