#pragma once

#include "neurox/neurox.h"

/**
 * @brief The Circuit class
 * Represents all biological data on the Global Address Space
 */
class Brain
{
  public:

    Brain()=delete;
    Brain(const int neuronsCount,
          const hpx_t neuronsAddr, const Mechanism * mechsTypes,
          const size_t mechsTypesCount, const int * mechDependencies);
    ~Brain();

    //global vars: all localities hold the same value
    int neuronsCount; 	///> total neurons count in the system
    hpx_t neuronsAddr; 	///> hpx address of the first position of the neurons array

    //Mechanisms data:
    int mechsTypesCount;    ///> number of mechanisms
    Mechanism * mechsTypes; ///> Unique information per mechanism type

    //BeforeAfter function list (few functions per mechanism)
    mod_f_t *beforeAfterFunction[BEFORE_AFTER_SIZE];

    static void registerHpxActions(); ///> Registers all HPX actions
    static hpx_action_t init; ///> Initializes Circuit as a copy
    static hpx_action_t clear; ///> deletes all data (including neurons and branches)
    static hpx_action_t finitialize; ///> finitialize.c::nrn_finitialize
    static hpx_action_t solve; ///>netpar.cpp:BBS_netpar_solve()

    inline hpx_t getNeuronAddr(int i) const {
        return hpx_addr_add(neuronsAddr, sizeof(Neuron)*i, sizeof(Neuron));
    }; ///> returns hpx address for i-th neuron
  private:

    static int callBeforeAfterMethod(const int functionId);
    static int finitialize_handler(); ///>finitialize.c::finitialize()
    static int solve_handler(); ///>netpar.cpp:BBS_netpar_solve()
    static int clear_handler(); ///> HPX destructor
    static int init_handler(const int neuronsCount, const hpx_t neuronsAddr,
                            const Mechanism *mechsTypes, const size_t mechsTypesCount,
                            const int *mechDependencies) ; ///>HPX constructor
} ;
