#pragma once

#include "neurox/neurox.h"
#include <array>

namespace neurox
{

/**
 * @brief The VecPlayContinuouX class
 * Equivalent to Coreneuron's VecPlayContinuous in nrniv/vrecitem.h
 * Includes/stores (?) information about continuous events
 */
class VecPlayContinuousX : Event
{
  public:
    VecPlayContinuousX() = delete;
    VecPlayContinuousX(double*, IvocVect* yvec, IvocVect* tvec, IvocVect* discon);
    ~VecPlayContinuousX();

    void deliver(floble_t t, Branch* branch); ///> at associated DiscreteEvent
    void continuous(floble_t tt); ///> play - every f(y, t) or res(y', y, t); record - advance_tn and initialize flag
    void play_init(Branch * branch); ///> called near beginning of finitialize

    //methods from VectPlayContinuos (nrniv/vrecitem.h)
    double interpolate(double tt);
    double interp(double th, double x0, double x1) {
        return x0 + (x1 - x0) * th;
    }
    void search(double tt);

    virtual int type() {
        return VecPlayContinuousType;
    }

    //PlayRecord vars:
    floble_t *pd_; ///>index to last value pointed by function continuous()

    //VecPlayContinuoys vars
    IvocVect* y_;
    IvocVect* t_;
    IvocVect* discon_indices_;
    size_t last_index_;
    size_t discon_index_;
    size_t ubound_index_;
};

}
