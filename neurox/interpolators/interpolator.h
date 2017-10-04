#pragma once

#include "neurox.h"

namespace neurox {

namespace interpolators {

/**
 * @brief The Interpolators enum
 * Identifies the algorithm and jacobian used on
 * fixed or variable step interpolations;
 */
enum class Interpolators: int {
  kCvodesNeuronSolver = 0,
  kCvodesDenseMatrix = 1,
  kCvodesDiagMatrix = 2,
  kCvodesSparseMatrix = 3,
  kBackwardEuler = 9
};

};  // interpolators

};  // neurox

// TODO can we move this somewhere else?
#include "neurox/interpolators/variable_time_step.h"
//#include "neurox/interpolators/backward_euler.h"
