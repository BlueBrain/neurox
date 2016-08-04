#include "neurox/Neurox.h"
#include <numeric>
#include <algorithm>

using namespace Neurox;
using namespace Neurox::Solver;

BackwardEuler::~BackwardEuler() {}

void BackwardEuler::solve(double dt, double tstop)
{
    //for all neurons
    hpx_par_for_sync( [&] (int i, void*) -> int
        {
            double dt = inputParams->dt;
            double commStepSize = 0.1; //TODO get the collective minimum value

            //loop through every step
            for (double t=0; t<inputParams->tstop; t += commStepSize )
            {
                //computation steps
                for (double tc=0; tc<commStepSize; tc += dt)
                    hpx_call_sync(getNeuronAddr(i), BackwardEuler::step, NULL, 0);

                //communication step (spike exchange)
                //TODO [...]
            }

            return HPX_SUCCESS;
        },
    0, neuronsCount, NULL);
}

hpx_action_t BackwardEuler::step = 0;
int BackwardEuler::step_handler()
{
    neurox_hpx_pin(Neuron); //We are in a Neuron

    //1. multicore.c::nrn_thread_table_check()
    local->callModFunction(Mechanism::ModFunction::threadTableCheck);

    //2. multicore.c::nrn_fixed_step_thread()

    //2.1 cvodestb::deliver_net_events(nth);
    //(send outgoing spikes netcvode.cpp::NetCvode::check_thresh() )
    //TODO check APT: if (local->topBranch->v[0] local->thresholdAP)
    bool reachedThresold = local->getSomaVoltage() >= local->APthreshold;
    hpx_addr_t spikesLco = reachedThresold ? Neuron::fireActionPotential(local) : HPX_NULL;

    //3. netcvode.cpp::NetCon::deliver()
    //calls NET_RECEIVE in mod files to receive synapses
    //netcvode.cpp::NetCvode::deliver_net_events()
    //            ->NetCvode::deliver_events()
    //            ->NetCvode::deliver_event()
    //            ->NetCon::deliver()
    //            ->net_receive() on mod files
    local->callNetReceiveFunction(0);

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
    hpx_call_sync(local->soma, Branch::callModFunction, NULL, 0, Mechanism::ModFunction::capacityCurrent, local->t, local->dt);
    //TODO: this is not a MOD file function, its in capac.c, has to be converted!

    local->t += .5*dt;

    //	fixed_play_continuous(nth); TODO
    hpx_call_sync(local->soma, Branch::callModFunction, NULL, 0, Mechanism::ModFunction::state, local->t, local->dt); //nonvint(nth);
    hpx_call_sync(local->soma, Branch::callModFunction, NULL, 0, Mechanism::ModFunction::after_solve, local->t, local->dt);

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
    neurox_hpx_register_action(0, step);
}
