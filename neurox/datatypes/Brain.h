#pragma once

#include "neurox/neurox.h"

/**
 * @brief The Circuit class
 * Represents all biological data on the Global Address Space
 */
class Brain
{
  public:

    Brain();
    ~Brain();

    //global vars: all localities hold the same value
    int neuronsCount; 	///> total neurons count in the system
    hpx_t neuronsAddr; 	///> hpx address of the first position of the neurons array

    //Mechanisms data:
    int mechanismsCount;    ///> number of mechanisms
    Mechanism * mechanisms; ///> Unique information per mechanism type

    static void registerHpxActions();		///> Registers all HPX actions
    static hpx_action_t initialize; ///> Initializes Circuit as a copy
    static hpx_action_t clear; ///> deletes all data (including neurons and branches)

    inline hpx_t operator [](int i) const {
        return hpx_addr_add(neuronsAddr, sizeof(Neuron)*i, sizeof(Neuron));
    };

  private:
    static int initialize_handler(const Brain * brain_new,
                                  const Mechanism * mechanisms, const int mechanismsCount,
                                  const int * mechDependencies, const int totalDependenciesCount);

    static int clear_handler();
} ;
