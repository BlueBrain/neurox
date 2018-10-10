#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::interpolators;

Interpolator* Interpolator::New(InterpolatorIds type, void * mem_addr) {
  switch (type) {
    case InterpolatorIds::kCvodeDenseMatrix:
      return mem_addr ? new(mem_addr) VariableTimeStep() : new VariableTimeStep();
    case InterpolatorIds::kCvodeDiagonalMatrix:
      return mem_addr ? new(mem_addr) VariableTimeStep() : new VariableTimeStep();
    case InterpolatorIds::kCvodePreConditionedDiagSolver:
      return mem_addr ? new(mem_addr) VariableTimeStep() : new VariableTimeStep();
    case InterpolatorIds::kCvodeSparseMatrix:
      return mem_addr ? new(mem_addr) VariableTimeStep() : new VariableTimeStep();
    case InterpolatorIds::kBackwardEuler:
      return mem_addr ? new(mem_addr) BackwardEuler() : new BackwardEuler();
    default:
      return nullptr;
  }
  assert(0);
  return nullptr;
}

size_t Interpolator::Size(InterpolatorIds type) {
  switch (type) {
    case InterpolatorIds::kCvodeDenseMatrix:
      return sizeof(VariableTimeStep);
    case InterpolatorIds::kCvodeDiagonalMatrix:
      return sizeof(VariableTimeStep);
    case InterpolatorIds::kCvodePreConditionedDiagSolver:
      return sizeof(VariableTimeStep);
    case InterpolatorIds::kCvodeSparseMatrix:
      return sizeof(VariableTimeStep);
    case InterpolatorIds::kBackwardEuler:
      return sizeof(BackwardEuler);
    default:
      assert(0);
      return -1;
  }
  assert(0);
  return -1;
}
