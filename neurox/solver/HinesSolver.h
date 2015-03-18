#pragma once

#include <map>
#include <queue>
#include <vector>
#include "neurox/neurox.h"

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

  static void SynchronizeThresholdV(Branch *local, floble_t *thresholdV = NULL);
  static void ResetMatrixRHSandD(Branch *local);
  static void SetupMatrixRHS(Branch *local);
  static void SetupMatrixDiagonal(Branch *local);
  static void BackwardTriangulation(Branch *local);
  static void ForwardSubstituion(Branch *local);
  static void UpdateV(Branch *local);

 private:
};

};  // namespace
};  // namespace
