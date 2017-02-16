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

    static void resetMatrixRHSandD(Branch * local);
    static void setupMatrixRHS(Branch * local);
    static void setupMatrixDiagonal(Branch * local);
    static void backwardTriangulation(Branch *local);
    static void forwardSubstituion(Branch *local);

  private:
};

}; //namespace
}; //namespace
