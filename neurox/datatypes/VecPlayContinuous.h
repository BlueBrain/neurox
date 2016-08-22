#pragma once

#include "neurox/Neurox.h"
#include <vector>

namespace NeuroX
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
    VecPlayContinuouX(double *pd, double *t, double *y, size_t size);
    ~VecPlayContinuouX();

    void deliver(double t, Branch* branch) override;
    void continuous(double tt);
    double getFirstInstant();

    double *pd; ///>index to last value pointed by function continuous()

  private:
    double * y; //array of values
    double * t; //array of times instantes
    int size;
    int uBoundIndex;
    int lastIndex;
};

}
