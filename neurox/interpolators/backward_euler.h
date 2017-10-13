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
class BackwardEuler : public Interpolator {

 public:

  BackwardEuler() {}

  const char* GetString() override;
  void Init(Branch*) override;
  void StepTo(Branch*, const double tend) override;

  static hpx_action_t Finitialize;

  /// HPX actions registration
  static void RegisterHpxActions();

  static void Finitialize2(Branch*); ///> finitialize.c::finitialize()
  static void Step(Branch *);

  static void SetupTreeMatrix(Branch*);
  static int GetTotalStepsCount();

 private:

  static int Finitialize_handler();

};

};  // namespace
};  // namespace
