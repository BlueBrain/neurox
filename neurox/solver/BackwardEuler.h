#pragma once

#include "neurox/Neurox.h"

namespace NeuroX
{

namespace Solver
{
/**
 * @brief The BackwardEuler class
 * Includes all logic for Bulk Synchronous Parallel execution of Backward Euler
 */
class BackwardEuler
{
  public:
    BackwardEuler()=delete;
    ~BackwardEuler();

    static void registerHpxActions(); ///> Registers all HPX actions
    static void run(); ///> netpar.cpp:BBS_netpar_solve()
    static hpx_action_t initialize;  ///> finitialize.c
    static hpx_action_t initializeBranch;  ///> finitialize.c
    static hpx_action_t solve; ///> performs one Backward Euler step in all neurons
    static hpx_action_t branchStep; ///> performs one Backward Euler step on the top branch

  private:
    static int solve_handler();
    static int branchStep_handler(const int * stepsCount_ptr, const size_t);
    static int initialize_handler(); ///> initialize.c
    static int initializeBranch_handler(const double * v, const size_t v_size);
} ;

} //Solver

} //NeuroX
