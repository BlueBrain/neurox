#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::interpolators;

Interpolator* Interpolator::New(InterpolatorIds type) {
  switch (type) {
    case InterpolatorIds::kCvodeDenseMatrix:
      return new VariableTimeStep();
    case InterpolatorIds::kCvodeDiagonalMatrix:
      return new VariableTimeStep();
    case InterpolatorIds::kCvodePreConditionedDiagSolver:
      return new VariableTimeStep();
    case InterpolatorIds::kCvodeSparseMatrix:
      return new VariableTimeStep();
    case InterpolatorIds::kBackwardEuler:
      return new BackwardEuler();
    default:
      return nullptr;
  }
  assert(0);
  return nullptr;
};
