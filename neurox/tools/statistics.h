/*
# =============================================================================
# Copyright (c) 2015 - 2021-2022 Blue Brain Project/EPFL
#
# See top-level LICENSE file for details.
# =============================================================================.
*/
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
  /// count of communication when cmd line flag "output comm counts" is enabled
  class CommCount {
   public:
    static struct Counts {
      unsigned point_to_point_count;
      unsigned reduce_count;
      unsigned spike_count;
    } counts;
    static hpx_t mutex;
    static hpx_t allreduce_lco;
    static hpx_t allreduce_future;
    static int allreduce_id;

    static hpx_action_t Init;
    static hpx_action_t Reduce;
    static hpx_action_t Subscribe;
    static hpx_action_t Join;
    static hpx_action_t Unsubscribe;
    static void Init_handler(Counts*, const size_t);
    static void Reduce_handler(Counts*, const Counts*, const size_t);
    static int Subscribe_handler(const hpx_t*, const size_t);
    static int Join_handler();
    static int Unsubscribe_handler(const hpx_t*, const size_t);

    static void ReduceCounts(Counts* s);
  };

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
