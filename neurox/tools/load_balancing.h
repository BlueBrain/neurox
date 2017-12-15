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
  static double GetWorkPerBranchSubsection(const double neuron_time,
                                           const int my_neurons_count);

  /// returns max work per locality, for branch parallelism
  static double GetWorkPerLocality(const double neuron_time,
                                   const int my_neurons_count);

  /// sums given runtime to variable total_mech_instances_runtime_
  static void AddToTotalMechInstancesRuntime(double);

  /// returns the workload assigned to each cluster of mech instances
  static double GetWorkloadPerMechInstancesThread();

  /// branch parallelism multiplier to compute max workload per section
  static constexpr double kBranchSectionsPerComputeUnit = 2;

  /// number of instances in cluster, for mechanism type parallelism
  static constexpr double kMechInstancesPercentagePerComputeUnit = 0.1;

  static void RegisterHpxActions();

 private:
  /// computation (ms) per compute node
  static double *load_balancing_table_;

  /// mutex for loadBalancingTable
  static hpx_t load_balancing_mutex_;

  /// total time spent on current+state mechanisms execution
  static double total_mech_instances_runtime_;
};

};  // Tools
};  // NeuroX
