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
    VecPlayContinuouX(floble_t *pd, const floble_t *t,
                      const floble_t *y, const size_t size);
    ~VecPlayContinuouX();

    void deliver(floble_t t, Branch* branch) override;
    void continuous(floble_t tt);
    floble_t getFirstInstant();
    int type() override { return VecPlayContinuousType; }

    floble_t *pd; ///>index to last value pointed by function continuous()

  private:
    floble_t * y; ///> array of values
    floble_t * t; ///> array of times instantes
    size_t size; ///> size of t and y arrays
    size_t uBoundIndex;
    size_t lastIndex;
};

}
