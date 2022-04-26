/*
# =============================================================================
# Copyright (c) 2015 - 2021-2022 Blue Brain Project/EPFL
#
# See top-level LICENSE file for details.
# =============================================================================.
*/
#pragma once

#include "neurox/neurox.h"

#include <map>
#include <queue>
#include <vector>

namespace neurox {

/**
 * @brief The Hines Solver class
 * Handles the Gaussian Elimination synchronizer according to Hines;
 */
class HinesSolver {
 public:
  HinesSolver() = delete;
  ~HinesSolver();

  /// Communicates A and B elements between branch subsections
  static void InitializeBranchConstants(const Branch *branch);

  /// Communicates voltage between soma and Axon Initial Segment
  static double GetAxonInitialSegmentVoltage(const Branch *branch);

  /// Zeros an array @arr for the given branch @branch
  static void ResetArray(const Branch *branch, floble_t *arr);

  /// Set-up of Matrix Right-Hand Side (pre-Gaussian)
  static void SetupMatrixRHS(Branch *branch);

  /// Set-up of Matrix Diagonal (pre-Gaussian)
  static void SetupMatrixDiagonal(Branch *branch);

  /// The Backward Triangulation step in inverted Gaussian Elimination
  static void BackwardTriangulation(Branch *branch);

  /// The Forward Substitution step in inverted Gaussian Elimination
  static void ForwardSubstitution(Branch *branch);

  /// Calls Triangulation and Substitution methods
  static void SolveTreeMatrix(Branch *branch);

  /// Update Voltage contributions from the Right-Hand Side
  static void UpdateVoltagesWithRHS(Branch *branch);

  /// Update Voltage contributions from branching points
  static void UpdateBranchVoltagesWithRHS(Branch *branch);

  // CVODE-specific methods
  static void ResetRHSandDNoCapacitors(Branch *);
  static void ResetRHSNoCapacitors(Branch *);
  static void SetupMatrixVoltageNoCapacitors(Branch *);

 private:
};

};  // namespace neurox
