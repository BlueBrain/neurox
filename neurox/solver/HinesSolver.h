#pragma once

#include "neurox/Neurox.h"
#include <queue>
#include <map>
#include <vector>

namespace  NeuroX {


namespace Solver {

/**
 * @brief The Hines Solver class
 * Handles the Gaussian Elimination algorithm according to Hines;
 */
class HinesSolver
{
  public:
    HinesSolver()=delete; ///> no constructor, build using hpx init function instead
    ~HinesSolver();

    static void registerHpxActions(); ///> Register all HPX actions
    static hpx_action_t setupMatrixRHS; ///> finitialize.c::nrn_finitialize
    static hpx_action_t setupMatrixLHS; ///> finitialize.c::nrn_finitialize
    static hpx_action_t gaussianBackTriangulation; ///> Gaussian elimination's back triangulation: solve_core.c:triang()
    static hpx_action_t gaussianFwdSubstitution; ///> Gaussian elimination's forward substitution: solve_core.c:bksub()

  private:

    static int setupMatrixRHS_handler(const double * parentV_ptr, const size_t);
    static int setupMatrixLHS_handler();
    static int gaussianBackTriangulation_handler();
    static int gaussianFwdSubstitution_handler(const double *, const size_t);
};

}; //namespace
}; //namespace
