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
  class SizeInfo {
   public:
    SizeInfo();
    ~SizeInfo();

    neuron_id_t neuron_id_;
    double morphologies_;
    double mechanisms_;
    double synapses_;
    double metadata_;
    double global_vars_;
    unsigned long long compartments_count_;
    unsigned long long branches_count_;
    unsigned long long mechs_instances_count_;

    double getTotalSize();

    SizeInfo& operator+=(const SizeInfo& rhs);
  };

  /// Register all HPX actions
  static void RegisterHpxActions();

  /// returns total simulation size
  static void OutputSimulationSize(bool write_to_file = true);

  /// print the statistics about mechanisms distribution per type
  static void OutputMechanismsDistribution(bool write_to_file = true);

  /// returns branch size, including branches, in KB
  static hpx_action_t GetNeuronSize;

  /// returns mechanisms count per type
  static hpx_action_t GetNeuronMechanismsDistribution;

 private:
  static int GetNeuronSize_handler();
  static int GetNeuronMechanismsDistribution_handler();
};
}  // namespace tools
}  // namespace neurox
