#include "neurox/neurox.h"
#include <numeric>
#include <algorithm>
#include <math.h>

using namespace neurox;
using namespace neurox::Solver;

BackwardEuler::~BackwardEuler() {}

void BackwardEuler::run()
{
  hpx_par_for_sync( [&] (int i, void*) -> int //for all neurons
  {
     //finitialize.c (nrn_finitialize( 1, inputParams.voltage )
     if (HPX_LOCALITY_ID ==0 && i==0 ) printf("Neuron::initialize...\n");
     return hpx_call_sync(getNeuronAddr(i), BackwardEuler::initializeNeuron, NULL, 0);
     if (HPX_LOCALITY_ID ==0 && i==0) printf("Neuron::solve...\n");
     return hpx_call_sync(getNeuronAddr(i), BackwardEuler::solve, NULL, 0);
     return HPX_SUCCESS;
  },
  0, neuronsCount, NULL);
}

hpx_action_t BackwardEuler::initializeBranch = 0;
int BackwardEuler::initializeBranch_handler(const double * v, const size_t v_size)
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(BackwardEuler::initializeBranch, v, v_size);
    //nrn_spike_exchange_init() inserts NetParEvent on queue;
    //Not applicable for our use case
    //nrn_play_init() initiates the VectorPlayContinuous
    local->initEventsQueue();
    local->deliverEvents_handler(NULL, 0);

    //set up by finitialize.c:nrn_finitialize(): if (setv)
    for (int n=0; n<local->n; n++)
        local->v[n]=*v;
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t BackwardEuler::initializeNeuron=0;
int BackwardEuler::initializeNeuron_handler()
{
    neurox_hpx_pin(Neuron);

    //set up by finitialize.c:nrn_finitialize(): if (setv)
    double v = inputParams->voltage;
    hpx_call_sync(local->soma, BackwardEuler::initializeBranch, NULL, 0, &v, sizeof(v));

    // the INITIAL blocks are ordered so that mechanisms that write
    // concentrations are after ions and before mechanisms that read
    // concentrations.
    local->callModFunction(Mechanism::ModFunction::before_initialize);
    local->callModFunction(Mechanism::ModFunction::initialize);
    local->callModFunction(Mechanism::ModFunction::after_initialize);

    //set up by finitialize.c:nrn_finitialize() -> fadvance_core.c:dt2thread()
    //local->cj = inputParams->secondorder ? 2.0/inputParams->dt : 1.0/inputParams->dt;
    //done when calling mechanisms //TODO have a copy per branch to speed-up?

    hpx_call_sync(local->soma, Branch::deliverEvents, NULL, 0, NULL,0);
    local->setupTreeMatrixMinimal();
    hpx_call_sync(local->soma, Branch::deliverEvents, NULL, 0, NULL, 0);

    //part of nrn_fixed_step_group_minimal
    //1. multicore.c::nrn_thread_table_check()
    local->callModFunction(Mechanism::ModFunction::threadTableCheck);

    neurox_hpx_unpin;
}


hpx_action_t BackwardEuler::branchStep = 0;
int BackwardEuler::branchStep_handler(const int  * stepsCount_ptr, const size_t size)
{
    /*
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(BackwardEuler::branchStep, stepsCount_ptr, size);

    for (int s=0; s<*stepsCount_ptr; s++)
    {
        //1. multicore.c::nrn_thread_table_check()
        Mechanism::ModFunction func;
        func = Mechanism::ModFunction::threadTableCheck;
        callModFunction_handler(&func, sizeof(func));

        //2.1 cvodestb::deliver_net_events(nth);
        //(send outgoing spikes netcvode.cpp::NetCvode::check_thresh() )
        bool reachedThresold = local->isSoma() && local->v[0]>APthreshold;
        hpx_addr_t spikesLco = reachedThresold ? local->fireActionPotential() : HPX_NULL;

    }

    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
    */
}

hpx_action_t BackwardEuler::solve = 0;
int BackwardEuler::solve_handler()
{
  neurox_hpx_pin(Neuron); //We are in a Neuron

  double dt = inputParams->dt;
  double dt_comm = 0.1; //TODO get the collective minimum value
  int stepsCountPerComm = dt_comm / dt;

  //this is void nrn_fixed_step_group_minimal(int n) in fadvance_core.c
  for (local->t=0; local->t<inputParams->tstop; local->t+=dt)
  {
    //2. multicore.c::nrn_fixed_step_thread()
    //2.1 cvodestb::deliver_net_events(nth);
    //(send outgoing spikes netcvode.cpp::NetCvode::check_thresh() )
    bool reachedThresold = local->getSomaVoltage() >= local->APthreshold;
    if (reachedThresold) local->fireActionPotential();
    local->t += .5*dt;
    hpx_call_sync(local->soma, Branch::deliverEvents, NULL, 0, &local->t, sizeof(local->t));
    hpx_call_sync(local->soma, Branch::fixedPlayContinuous, NULL, 0);

    local->setupTreeMatrixMinimal();

    //eion.c : second_order_cur()
    if (inputParams->secondorder == 2)
        hpx_call_sync(local->soma, Branch::secondOrderCurrent, NULL, 0);

    //fadvance_core.c : update()
    hpx_call_sync(local->soma, Branch::updateV, NULL, 0, inputParams->secondorder, sizeof(inputParams->secondorder));
    local->callModFunction(Mechanism::ModFunction::jacob);
    local->t += .5*dt;
    hpx_call_sync(local->soma, Branch::fixedPlayContinuous, NULL, 0);
    local->callModFunction(Mechanism::ModFunction::state); //nonvint
    local->callModFunction(Mechanism::ModFunction::after_solve);
    hpx_call_sync(local->soma, Branch::deliverEvents, NULL, 0);

    //if we are at the output time instant output to file
    if (fmod(local->t, inputParams->dt_io) == 0)
    {
        //output
    }

    //make sure all synapses from N steps before were delivered
    //(thus other neurons wait for this neuron one before stepping)
    local->waitForSynapsesDelivery(stepsCountPerComm);

    //make sure I'm not more than N steps ahead of that connect to me!
    //TODO
  }
  neurox_hpx_unpin;
}

void BackwardEuler::registerHpxActions()
{
    neurox_hpx_register_action(0, BackwardEuler::solve);
    neurox_hpx_register_action(1, BackwardEuler::branchStep);
    neurox_hpx_register_action(0, BackwardEuler::initializeNeuron);
    neurox_hpx_register_action(1, BackwardEuler::initializeBranch)
}
