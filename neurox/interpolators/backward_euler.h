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
  hpx_t StepTo(Branch* branch, const double tstop) override;

  static hpx_t Step(Branch*, const bool benchmark = false);

  static hpx_t StepParallel(Branch*);

  static hpx_action_t Finitialize;

  static void Finitialize2(Branch*);  ///> finitialize.c::finitialize()
  static void SetupTreeMatrix(Branch*);
  static int GetTotalSteps();
  static int GetMinSynapticDelaySteps();

  inline static void RunFunction2(hpx_t, const Branch*, const int);

  static hpx_action_t RunFunction;
  static void RegisterHpxActions();
 private:
  static int RunFunction_handler(const int, const void* [], const size_t[]);
};

};  // namespace
};  // namespace
