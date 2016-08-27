#pragma once

#include "neurox/neurox.h"
#include <queue>
#include <map>
#include <vector>

namespace  neurox {


namespace Solver {

/**
 * @brief The Hines Solver class
 * Handles the Gaussian Elimination algorithm according to Hines;
 */
class HinesSolver
{
  public:
    HinesSolver()=delete;
    ~HinesSolver();

    static void gaussianFwdTriangulation(Branch * local);
    static void gaussianBackSubstitution(Branch * local);
  private:
};

}; //namespace
}; //namespace
