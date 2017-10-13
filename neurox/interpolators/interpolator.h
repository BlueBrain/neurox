#pragma once

#include "neurox.h"

namespace neurox {

namespace interpolators {

/**
 * @brief The Interpolators enum
 * Identifies the synchronizer and jacobian used on
 * fixed or variable step interpolations;
 */
enum class Interpolators : int {
  kCvodePreConditionedDiagSolver = 0,
  kCvodeDenseMatrix = 1,
  kCvodeDiagonalMatrix = 2,
  kCvodeSparseMatrix = 3,
  kBackwardEuler = 9
};

class Interpolator {

 public:

    /// Returns class type as string
    const virtual char* GetString() = 0;

    virtual void Init(Branch*) {}
    virtual void StepTo(Branch*, const double tend) = 0;
    virtual void Clear(Branch*) {}

    /// Returns an instantiated class of the given type
    static Interpolator* New(Interpolators);

};
};  // interpolators

};  // neurox
