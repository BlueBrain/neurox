#pragma once

#include "neurox/neurox.h"

namespace neurox{
namespace tools{

/**
 * @brief The Neuron class
 * Represents a neuron as its metadata and a tree-based orphology
 */
class Statistics
{
  public:
    static void RegisterHpxActions(); ///> Register all HPX actions
    static void OutputSimulationSize(bool writeToFile = true); ///> returns total simulation size
    static void OutputMechanismsDistribution(bool writeToFile = true); ///> print the statistics about mechanisms distribution per type
    static hpx_action_t GetNeuronSize; ///> returns branch size, including branches, in KB
    static hpx_action_t GetNeuronMechanismsDistribution; ///> returns mechanisms count per type

  private:
    static int GetNeuronSize_handler();
    static int GetNeuronMechanismsDistribution_handler();
    class SizeInfo;
};
} //Tools
} //neurox
