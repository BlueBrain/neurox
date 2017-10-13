#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::interpolators;

Interpolator* Interpolator::New(Interpolators type) {
  switch (type) {
    case Interpolators::kCvodeDenseMatrix:
      return new VariableTimeStep();
    case Interpolators::kCvodeDiagonalMatrix:
      return new VariableTimeStep();
    case Interpolators::kCvodePreConditionedDiagSolver:
      return new VariableTimeStep();
    case Interpolators::kCvodeSparseMatrix:
      return new VariableTimeStep();
    case Interpolators::kBackwardEuler:
      return new BackwardEuler();
    default:
      return nullptr;
  }
  assert(0);
  return nullptr;
};
