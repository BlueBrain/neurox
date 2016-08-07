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
    static void printSimulationSize();  ///> returns total simulation size
    static hpx_action_t printNeuronSize;      ///> returns neuron size, including branches, in KB
    static hpx_action_t printBranchSize;      ///> returns branch size, including branches, in KB

  private:
    static int printNeuronSize_handler();
    static int printBranchSize_handler();
    class SizeInfo;
};
}
}
