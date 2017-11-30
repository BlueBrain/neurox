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

  /// Queries/updates load balancing table, for load balancing
  static int QueryLoadBalancingTable_handler(const int nargs,
                                             const void *args[],
                                             const size_t[]);
  /// Prints load balancing tabl
  static void PrintLoadBalancingTable();

  /// returns work per branch-subsection, for branch parallelism
  static double GetWorkPerBranchSubsection(const double neuron_time,
                                           const int neurons_count);

  /// branch parallelism multiplier to compute max workload per section
  static constexpr double kSubSectionsPerComputeUnit = 2;

  /// number of instances in cluster, for mechanism type parallelism
  static const int kMechInstancesPerCluster = 1000;

  static void RegisterHpxActions();

 private:
  /// computation (ms) per compute node
  static double *load_balancing_table_;

  /// mutex for loadBalancingTable
  static hpx_t load_balancing_mutex_;
};

};  // Tools
};  // NeuroX
