/*
# =============================================================================
# Copyright (c) 2015 - 2021-2022 Blue Brain Project/EPFL
#
# See top-level LICENSE file for details.
# =============================================================================.
*/
#pragma once

#include "neurox/neurox.h"

#include <deque>
#include <map>
#include <memory>
#include <vector>

using namespace std;

namespace neurox {

namespace tools {

/**
 * @brief The LoadBalancing class
 */
class LoadBalancing {
 public:
  LoadBalancing();
  ~LoadBalancing();

  static hpx_action_t QueryLoadBalancingTable;
  static hpx_action_t UpdateLoadBalancingTable;

  /// Queries load balancing table for least busy compute node
  static int QueryLoadBalancingTable_handler();

  /// Updates load balancing table with locality and runtime
  static int UpdateLoadBalancingTable_handler(const int nargs,
                                              const void *args[],
                                              const size_t[]);

  /// Prints load balancing tabl
  static void PrintLoadBalancingTable();

  /// returns work per branch-subsection, for branch parallelism
  static double GetMaxWorkPerBranchSubTree(const double neuron_time,
                                           const int my_neurons_count);

  /// returns max work per locality, for branch parallelism
  static double GetMaxWorkPerBranchSubSection(const double neuron_time,
                                              const int my_neurons_count);

  /// returns the workload assigned to each cluster of mech instances
  static double GetWorkloadPerMechInstancesBlock(
      const double total_mech_instances_runtime);

  /// register all HPX actions related to this class;
  static void RegisterHpxActions();

 private:
  /// computation (ms) per compute node
  static double *load_balancing_table_;

  /// mutex for loadBalancingTable
  static hpx_t load_balancing_mutex_;
};

};  // namespace tools
};  // namespace neurox
