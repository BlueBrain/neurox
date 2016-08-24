#include "neurox/neurox.h"
#include <algorithm>
#include <cstring>

using namespace neurox;

VecPlayContinuouX::VecPlayContinuouX(double * const pd, const double *t,
                                     const double *y, const size_t size)
    :size(size), uBoundIndex(0), lastIndex(0), pd(pd)
{
    assert(size>0);
    this->t = new double[size];
    this->y = new double[size];
    std::memcpy(this->t, t, size*sizeof(double));
    std::memcpy(this->y, y, size*sizeof(double));
}

VecPlayContinuouX::~VecPlayContinuouX()
{
    delete [] t;
    delete [] y;
};

void VecPlayContinuouX::continuous(double tt)
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

        const double x0 = y[lastIndex-1];
        const double x1 = y[lastIndex];
        const double t0 = t[lastIndex - 1];
        const double t1 = t[lastIndex];

        /*** alternative
        auto iter = lower_bound(t.begin(), t.end(), tt);
        int index = std::distance(t.begin(), iter);
        double x0 = y[index];
        double x1 = y[index+1];
        double t0 = t[index];
        double t1 = t[index+1];
        ***/

        if (t0 == t1)
            *pd = (x0 + x1)/2.;
        else //interp(...)
            *pd = x0 + (x1 - x0)*(tt - t0)/(t1 - t0);
    }
}

double  VecPlayContinuouX::getFirstInstant()
{
    return this->t[0];
}

void VecPlayContinuouX::deliver(double t, Branch* branch)
{
    if (uBoundIndex < this->size-1)
    {
        double nextDeliveryTime = this->t[++uBoundIndex];
        hpx_lco_sema_p(branch->eventsQueueMutex);
        branch->eventsQueue.push(make_pair(nextDeliveryTime, (Event*) this));
        hpx_lco_sema_v_sync(branch->eventsQueueMutex);
    }
    continuous(t);
}

