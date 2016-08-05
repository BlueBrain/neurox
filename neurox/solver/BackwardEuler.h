#pragma once

#include "neurox/Neurox.h"

namespace NeuroX
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
    ~BackwardEuler();

    static void registerHpxActions(); ///> Registers all HPX actions
    static void solve(); ///> netpar.cpp:BBS_netpar_solve()
    static hpx_action_t step; ///> performs one Backward Euler step in all neurons

  private:
    static int step_handler(const int * stepsCount_ptr, const size_t);
} ;

} //Solver

} //NeuroX
