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
 * Handles the Gaussian Elimination synchronizer according to Hines;
 */
class HinesSolver {
 public:
  HinesSolver() = delete;
  ~HinesSolver();

  static void SynchronizeThresholdV(const Branch *branch,
                                    floble_t *threshold_v = NULL);
  static void ResetArray(const Branch *branch, floble_t *arr);
  static void SetupMatrixRHS(Branch *branch);
  static void SetupMatrixDiagonal(Branch *branch);
  static void BackwardTriangulation(Branch *branch);
  static void ForwardSubstituion(Branch *branch);
  static void SolveTreeMatrix(Branch *branch);
  static void UpdateVoltagesWithRHS(Branch *branch);

  // CVODE-specific methods
  static void ResetRHSandDNoCapacitors(Branch *, void *);
  static void ResetRHSNoCapacitors(Branch *, void *);
  static void SetupMatrixVoltageNoCapacitors(Branch *, void *);

 private:
};

};  // namespace
};  // namespace
