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
    int multiSplit; 	///> 0 or 1 for multisplit or not
    hpx_t neuronsAddr; 	///> hpx address of the first position of the neurons array

    Mechanism * mechanisms; ///>Unique information per mechanism

    static void registerHpxActions();		///> Registers all HPX actions
    static hpx_action_t initialize; ///> Initializes Circuit as a copy

    inline hpx_t operator [](int i) const {
        return hpx_addr_add(neuronsAddr, sizeof(Neuron)*i, sizeof(Neuron));
    };

  private:
    static int initialize_handler(const Brain * circuit, const size_t size);
} ;
