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
          const hpx_t neuronsAddr, const Mechanism * mechanisms,
          const size_t mechanismsCount, const int * mechDependencies);
    ~Brain();

    //global vars: all localities hold the same value
    int neuronsCount; 	///> total neurons count in the system
    hpx_t neuronsAddr; 	///> hpx address of the first position of the neurons array

    //Mechanisms data:
    int mechanismsCount;    ///> number of mechanisms
    Mechanism * mechanisms; ///> Unique information per mechanism type

    static void registerHpxActions(); ///> Registers all HPX actions
    static hpx_action_t initialize; ///> Initializes Circuit as a copy
    static hpx_action_t clear; ///> deletes all data (including neurons and branches)

    inline hpx_t getNeuronAddr(int i) const {
        return hpx_addr_add(neuronsAddr, sizeof(Neuron)*i, sizeof(Neuron));
    }; ///> returns hpx address for i-th neuron

  private:
    static int initialize_handler(const int neuronsCount,
                                  const hpx_t neuronsAddr, const Mechanism * mechanisms,
                                  const size_t mechanismsCount, const int * mechDependencies);
                                  ///>HPX constructor

    static int clear_handler(); ///> HPX destructor
} ;
