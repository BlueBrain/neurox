#pragma once

#include "neurox/neurox.h"

#include <functional>  //std::greater_equal
#include <map>
#include <new>  //placement new
#include <queue>
#include <vector>

using namespace neurox;

namespace neurox {

namespace interpolators {

/**
 * @brief The Branch class
 * Represents a branch as a continuous set of compartments (nodes)
 * and a bifurcation made of branches branchs spreading out from the last node;
 * Also represents mechanisms as a list of mechanisms and their applications
 * to several nodes of this branch.
 */
class BackwardEuler {
 public:

  static void SetupTreeMatrix(Branch*);
  static void Finitialize2(Branch*);
  static void BackwardEulerStep(Branch *);

  static hpx_action_t Finitialize;  ///> finitialize.c::finitialize()
  static hpx_action_t Run;
  static hpx_action_t RunOnLocality;

  static void RegisterHpxActions();  ///> Register all HPX actions

 private:
  static int Finitialize_handler();
  static int Run_handler(const int*, const size_t);
  static int RunOnLocality_handler(const int*, const size_t);
};

};  // namespace
};  // namespace
