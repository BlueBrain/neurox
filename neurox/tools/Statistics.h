#pragma once

#include "neurox/neurox.h"

namespace neurox{
namespace Tools{

/**
 * @brief The Neuron class
 * Represents a neuron as its metadata and a tree-based orphology
 */
class Statistics
{
  public:
    static void registerHpxActions(); ///> Register all HPX actions
    static void outputSimulationSize(bool writeToFile = true); ///> returns total simulation size
    static void outputMechanismsDistribution(bool writeToFile = true); ///> print the statistics about mechanisms distribution per type
    static hpx_action_t getNeuronSize; ///> returns branch size, including branches, in KB
    static hpx_action_t getNeuronMechanismsDistribution; ///> returns mechanisms count per type

  private:
    static int getNeuronSize_handler();
    static int getNeuronMechanismsDistribution_handler();
    class SizeInfo;
};
} //Tools
} //neurox
