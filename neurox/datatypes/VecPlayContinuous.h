#pragma once

#include "neurox/neurox.h"
#include <vector>

namespace neurox
{

/**
 * @brief The VecPlayContinuouX class
 * Equivalent to Coreneuron's VecPlayContinuous in nrniv/vrecitem.h
 * Includes/stores (?) information about continuous events
 */
class VecPlayContinuouX : Event
{
  public:
    VecPlayContinuouX() = delete;
    VecPlayContinuouX(double * const pd, const double *t, const double *y, const size_t size);
    ~VecPlayContinuouX();

    void deliver(double t, Branch* branch) override;
    void continuous(double tt);
    double getFirstInstant();
    int type() override { return VecPlayContinuousType; }

    double *pd; ///>index to last value pointed by function continuous()

  private:
    double * y; ///> array of values
    double * t; ///> array of times instantes
    int size; ///> size of t and y arrays
    int uBoundIndex;
    int lastIndex;
};

}
