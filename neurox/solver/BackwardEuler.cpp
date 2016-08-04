#include "neurox/Neurox.h"
#include <numeric>
#include <algorithm>

using namespace Neurox;
using namespace Neurox::Solver;

BackwardEuler::~BackwardEuler() {}

void BackwardEuler::solve()
{
    //for all neurons
    hpx_par_for_sync( [&] (int i, void*) -> int
        {
            double dt = inputParams->dt;
            double t_io = inputParams->dt_io;
            double dt_comm = 0.1; //TODO get the collective minimum value
            hpx_t neuronAddr= getNeuronAddr(i);

            //for every communication step
            for (double t_comm=0; t_comm<inputParams->tstop; t_comm += dt_comm )
            {
                printf("time = %.2f ms\n", t_comm);

                //communication step (spike exchange)
                //TODO [...]

                //for every computation step inside of the communication step
                for (double t=t_comm; t<t_comm+dt_comm; t += dt)
                {
                    hpx_call_sync(neuronAddr, BackwardEuler::step, NULL, 0);

                    //if we are at the output time instant (or passed it), output to file
                    if (t >= t_io)
                    {
                        //output
                        t_io += inputParams->dt_io;
                    }
                }
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
    //TODO this voltage access should be done on a local memory access, see conversation with Luke
    bool reachedThresold = local->getSomaVoltage() >= local->APthreshold;
    hpx_addr_t spikesLco = reachedThresold ? local->fireActionPotential() : HPX_NULL;

    //3. netcvode.cpp::NetCon::deliver()
    //calls NET_RECEIVE in mod files to receive synapses
    //netcvode.cpp::NetCvode::deliver_net_events()
    //            ->NetCvode::deliver_events()
    //            ->NetCvode::deliver_event()
    //            ->NetCon::deliver()
    //            ->net_receive() on mod files
    local->callNetReceiveFunction(0); //TODO replace 0/1 by Function::NetReceive and Function::NetReceiveInit

    local->t += .5*dt;

    //TODO: fixed_play_continuous; for PlayRecord (whole logic missing)

    local->setupTreeMatrixMinimal();

    //Linear Algebra: Gaussian elimination. solve_core.c:nrn_solve_minimal()
    double parentRHS=0;
    hpx_call_sync(local->soma, Branch::gaussianBackTriangulation, NULL, 0);
    hpx_call_sync(local->soma, Branch::gaussianFwdSubstitution,   NULL, 0, &parentRHS, sizeof(parentRHS));

    //eion.c : second_order_cur()
    if (inputParams->secondorder == 2)
        hpx_call_sync(local->soma, Branch::secondOrderCurrent, NULL, 0);

    //fadvance_core.c : update()
    hpx_call_sync(local->soma, Branch::updateV, NULL, 0, inputParams->secondorder, sizeof(inputParams->secondorder));
    local->callModFunction(Mechanism::ModFunction::capacityCurrent);
    //TODO: this is not a MOD file function, its in capac.c, has to be converted!

    local->t += .5*dt;

    //	fixed_play_continuous(nth); TODO
    local->callModFunction(Mechanism::ModFunction::state);
    local->callModFunction(Mechanism::ModFunction::after_solve);

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
