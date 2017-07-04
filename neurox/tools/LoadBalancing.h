#pragma once

#include "neurox/neurox.h"

#include <map>
#include <vector>
#include <memory>
#include <deque>

using namespace std;

namespace neurox {

namespace Tools {

/**
 * @brief The LoadBalancing class
 */
class LoadBalancing
{
  public:
    LoadBalancing();
    ~LoadBalancing();

    static hpx_action_t queryLoadBalancingTable;

    static void print();
    static void registerHpxActions();
    static int queryLoadBalancingTable_handler(const int nargs, const void *args[], const size_t[]);

  private:
    static double *loadBalancingTable; ///> computation (ms) per compute node
    static hpx_t loadBalancingMutex;;  ///> mutex for loadBalancingTable
};

}; //Misc
}; //NeuroX

