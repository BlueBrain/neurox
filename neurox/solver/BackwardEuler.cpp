#include "neurox/Neurox.h"
#include <numeric>
#include <algorithm>

using namespace Neurox;
using namespace Neurox::Solver;

BackwardEuler::~BackwardEuler() {}

void BackwardEuler::solve(double dt, double tstop)
{
    for (double t=0; t<tstop; t+=dt)
        hpx_par_for_sync( [&] (int i, void*) -> int
            { return hpx_call_sync(getNeuronAddr(i), BackwardEuler::step, NULL, 0); },
            0, neuronsCount, NULL);
}

hpx_action_t BackwardEuler::step = 0;
int BackwardEuler::step_handler()
{
    neurox_hpx_pin(Neuron); //We are in a Neuron

    //1. multicore.c::nrn_thread_table_check()
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::Functions::threadTableCheck, local->t, local->dt);

    //2. multicore.c::nrn_fixed_step_thread()

    //2.1 cvodestb::deliver_net_events(nth);
    //(send outgoing spikes netcvode.cpp::NetCvode::check_thresh() )
    //TODO check APT: if (local->topBranch->v[0] local->thresholdAP)
    bool reachedThresold = true;
    hpx_addr_t spikesLco = reachedThresold ? Neuron::fireActionPotential(local) : HPX_NULL;

    //3. netcvode.cpp::NetCon::deliver()
    //calls NET_RECEIVE in mod files to receive synapses
    //TODO: are these the same ones sent, or from the previous interval?
    //netcvode.cpp::NetCvode::deliver_net_events()
    //            ->NetCvode::deliver_events()
    //            ->NetCvode::deliver_event()
    //            ->NetCon::deliver()
    //            ->net_receive() on mod files
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::Functions::pntReceive, local->t, local->dt);

    local->t += .5*dt;

    //TODO: fixed_play_continuous; for PlayRecord (whole logic missing)

    Neuron::setupTreeMatrixMinimal(local);

    //Linear Algebra: Gaussian elimination. solve_core.c:nrn_solve_minimal()
    char isSoma=1;
    hpx_call_sync(local->soma, Branch::gaussianBackTriangulation, NULL, 0, isSoma);
    hpx_call_sync(local->soma, Branch::gaussianFwdSubstitution, NULL, 0, isSoma);

    //eion.c : second_order_cur()
    if (inputParams->secondorder == 2)
        hpx_call_sync(local->soma, Branch::secondOrderCurrent, NULL, 0);

    //fadvance_core.c : update()
    hpx_call_sync(local->soma, Branch::updateV, NULL, 0, inputParams->secondorder);
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::Functions::capacityCurrent, local->t, local->dt);
    //TODO: this is not a MOD file function, its in capac.c, has to be converted!

    local->t += .5*dt;

    //	fixed_play_continuous(nth); TODO
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::Functions::state, local->t, local->dt); //nonvint(nth);
    hpx_call_sync(local->soma, Branch::callMechsFunction, NULL, 0, Mechanism::Functions::after_solve, local->t, local->dt);

    //wait for all post-synaptic neurons to receive (not process) synapses
    if (spikesLco != HPX_NULL)
    {
        hpx_lco_wait(spikesLco);
        hpx_lco_delete(spikesLco, HPX_NULL);
    }
    neurox_hpx_unpin;
}

void BackwardEuler::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED, step, step_handler, HPX_DOUBLE);
}
