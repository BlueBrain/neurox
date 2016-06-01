#pragma once

#include "neurox/Neurox.h"

namespace Neurox
{

namespace Solver
{
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
    static void solve(); ///> netpar.cpp:BBS_netpar_solve()

    static hpx_action_t step;
    static hpx_action_t queueSpike;

  private:
    static int step_handler(const double dt); ///> fadvance_core.c::nrn_fixed_step_minimal()
    static int queueSpike_handler(const Synapse * syn, size_t size); ///> queueing of synapses
} ;

} //Solver

} //Neurox
