#include "neurox/neurox.h"
#include <algorithm>
#include <cstring>

using namespace neurox;

VecPlayContinuouX::VecPlayContinuouX(floble_t * const pd, const floble_t *t,
                                     const floble_t *y, const size_t size)
    :size(size), uBoundIndex(0), lastIndex(0), pd(pd)
{
    assert(size>0);
    this->t = new floble_t[size];
    this->y = new floble_t[size];
    for (size_t i=0; i<size; i++)
    {
        this->t[i] = (floble_t) t[i];
        this->y[i] = (floble_t) y[i];
    }

}

VecPlayContinuouX::~VecPlayContinuouX()
{
    delete [] t;
    delete [] y;
};

void VecPlayContinuouX::continuous(floble_t tt)
{
    //vrecord.cpp::VecPlayContinuous::interpolate()
    if (tt >= this->t[uBoundIndex])
    {
        lastIndex = uBoundIndex;
        if (uBoundIndex==0)
            *pd = this->y[0];
    }
    else if (tt <= this->t[0])
    {
        *pd = this->y[0];
    }
    else
    {
        //search
        while (tt <  this->t[lastIndex]) lastIndex--;
        while (tt >= this->t[lastIndex]) lastIndex++;

        const floble_t x0 = y[lastIndex-1];
        const floble_t x1 = y[lastIndex];
        const floble_t t0 = t[lastIndex - 1];
        const floble_t t1 = t[lastIndex];

        /*** alternative
        auto iter = lower_bound(t.begin(), t.end(), tt);
        int index = std::distance(t.begin(), iter);
        floble_t x0 = y[index];
        floble_t x1 = y[index+1];
        floble_t t0 = t[index];
        floble_t t1 = t[index+1];
        ***/

        if (t0 == t1)
            *pd = (x0 + x1)/2.;
        else //interp(...)
            *pd = x0 + (x1 - x0)*(tt - t0)/(t1 - t0);
    }
}

floble_t VecPlayContinuouX::getFirstInstant()
{
    return this->t[0];
}

void VecPlayContinuouX::deliver(floble_t t, Branch* branch)
{
    if (uBoundIndex < this->size-1)
    {
        floble_t nextDeliveryTime = this->t[++uBoundIndex];
        hpx_lco_sema_p(branch->eventsQueueMutex);
        branch->eventsQueue.push(make_pair(nextDeliveryTime, (Event*) this));
        hpx_lco_sema_v_sync(branch->eventsQueueMutex);
    }
    continuous(t);
}

