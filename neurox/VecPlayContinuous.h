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
    VecPlayContinuousX(double*, size_t size_, floble_t * yvec, floble_t * tvec, floble_t * discon);
    ~VecPlayContinuousX();

    void Deliver(floble_t t, Branch* branch); ///> at associated DiscreteEvent
    void Continuous(floble_t tt); ///> play - every f(y, t) or res(y', y, t); record - advance_tn and initialize flag
    void PlayInit(Branch * branch); ///> called near beginning of finitialize

    //methods from VectPlayContinuos (nrniv/vrecitem.h)
    double Interpolate(double tt);
    inline double Interp(double th, double x0, double x1) {
        return x0 + (x1 - x0) * th;
    }
    void Search(double tt);

    virtual int Type() {
        return VecPlayContinuousType;
    }

    //PlayRecord vars:
    floble_t *pd_; ///>index to last value pointed by function continuous()

    //VecPlayContinuoys vars
    floble_t  * y_;
    floble_t * t_;
    floble_t* discon_indices_;
    size_t size_; ///> size of arrays t_, y_
    size_t last_index_;
    size_t discon_index_;
    size_t ubound_index_;
};

}
