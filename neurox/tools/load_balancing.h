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

  static void PrintTable();
  static void RegisterHpxActions();
  static int QueryLoadBalancingTable_handler(const int nargs,
                                             const void *args[],
                                             const size_t[]);

 private:
  /// computation (ms) per compute node
  static double *loadBalancingTable;

  /// mutex for loadBalancingTable
  static hpx_t loadBalancingMutex;
};

};  // Tools
};  // NeuroX
