#pragma once

#include "neurox/Neurox.h"

namespace NeuroX{
namespace Misc{

/**
 * @brief The Neuron class
 * Represents a neuron as its metadata and a tree-based orphology
 */
class Statistics
{
  public:
    static void registerHpxActions(); ///> Register all HPX actions
    static void printSimulationSize(); ///> returns total simulation size
    static void printMechanismsDistribution(); ///> print the statistics about mechanisms distribution per type
    static hpx_action_t getNeuronSize; ///> returns neuron size, including branches, in KB
    static hpx_action_t getBranchSize; ///> returns branch size, including branches, in KB
    static hpx_action_t getNeuronMechanismsDistribution; ///> returns mechanisms count per type
    static hpx_action_t getBranchMechanismsDistribution; ///> returns mechanisms count per type

  private:
    static int getNeuronSize_handler();
    static int getBranchSize_handler();
    static int getNeuronMechanismsDistribution_handler();
    static int getBranchMechanismsDistribution_handler();
    class SizeInfo;
};
}
}
