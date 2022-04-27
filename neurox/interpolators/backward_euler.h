/*
# =============================================================================
# Copyright (c) 2015 - 2021-2022 Blue Brain Project/EPFL
#
# See top-level LICENSE file for details.
# =============================================================================.
*/
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
  void StepTo(Branch* branch, const floble_t tstop) override;

  static hpx_t Step(Branch*, const bool benchmark = false);

  static hpx_action_t Finitialize;

  static void Finitialize2(Branch*);  ///> finitialize.c::finitialize()
  static void SetupTreeMatrix(Branch*);

 private:
  static int Finitialize_handler();
};

};  // namespace interpolators
};  // namespace neurox
