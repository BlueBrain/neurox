#pragma once

#include "neurox/neurox.h"

namespace neurox {
namespace tools {

/**
 * @brief The Neuron class
 * Represents a neuron as its metadata and a tree-based orphology
 */
class Statistics {
 public:
  /// Register all HPX actions
  static void RegisterHpxActions();

  /// returns total simulation size
  static void OutputSimulationSize(bool writeToFile = true);

  /// print the statistics about mechanisms distribution per type
  static void OutputMechanismsDistribution(bool writeToFile = true);

  /// returns branch size, including branches, in KB
  static hpx_action_t GetNeuronSize;

  /// returns mechanisms count per type
  static hpx_action_t GetNeuronMechanismsDistribution;

 private:
  static int GetNeuronSize_handler();
  static int GetNeuronMechanismsDistribution_handler();
  class SizeInfo;
};
}  // Tools
}  // neurox
