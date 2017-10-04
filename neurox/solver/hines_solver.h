#pragma once

#include "neurox/neurox.h"

#include <map>
#include <queue>
#include <vector>

using namespace interpolators;

namespace neurox {

namespace solver {

/**
 * @brief The Hines Solver class
 * Handles the Gaussian Elimination algorithm according to Hines;
 */
class HinesSolver {
 public:
  HinesSolver() = delete;
  ~HinesSolver();

  static void SynchronizeThresholdV(Branch *local,
                                    floble_t *threshold_v = NULL);
  static void ResetMatrixRHSandD(Branch *local);
  static void ResetMatrixRHS(Branch *local);
  static void ResetMatrixV(Branch *local);
  static void SetupMatrixRHS(Branch *local);
  static void SetupMatrixDiagonal(Branch *local);
  static void BackwardTriangulation(Branch *local);
  static void ForwardSubstituion(Branch *local);
  static void UpdateVoltagesWithRHS(Branch *local);

  // CVODE-specific methods
  static void ResetRHSandDNoCapacitors(Branch*, void*);
  static void SetupMatrixRHSNoCapacitors(Branch*, void*);

 private:
};

};  // namespace
};  // namespace
