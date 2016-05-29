#pragma once

#include "neurox/neurox.h"

/**
 * @brief The Circuit class
 * Represents all biological data on the Global Address Space
 */
class BackwardEuler
{
  public:
    BackwardEuler()=delete;
    BackwardEuler(InputParams * inputParams);
    ~BackwardEuler();

    static void registerHpxActions(); ///> Registers all HPX actions
    static int solve(); ///> netpar.cpp:BBS_netpar_solve()

    double t;
    double dt;
    double tstop;

    static hpx_action_t solve;
    static hpx_action_t step;
    static hpx_action_t queueSpike;
    static hpx_action_t deliverSpike;

  private:
    static int step_handler(const double dt); ///> fadvance_core.c::nrn_fixed_step_minimal()
    static int queueSpike_handler(const Synapse * syn, size_t size); ///> queueing of synapses
    static int deliverSpike_handler(const int time); ///> delivering of (queued) synapses
} ;
