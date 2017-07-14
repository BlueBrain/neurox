#include "neurox/neurox.h"
#include <algorithm>
#include <cstring>

using namespace neurox;

VecPlayContinuousX::VecPlayContinuousX(double* pd, size_t size, floble_t * yvec, floble_t * tvec, floble_t * discon)
    :pd_(pd), size_(size), y_(yvec), t_(tvec), discon_indices_(discon), last_index_(0), discon_index_(0), ubound_index_(0)
{}

void VecPlayContinuousX::PlayInit(Branch * branch) {
    last_index_ = 0;
    discon_index_ = 0;
    if (discon_indices_) {
        if (size_ > 0) {
            ubound_index_ = (int)discon_indices_[discon_index_++];
            branch->AddEventToQueue( t_[ubound_index_], (Event*) this);
        } else {
            ubound_index_ = size_ - 1;
        }
    } else {
        ubound_index_ = 0;
        branch->AddEventToQueue(t_[ubound_index_], (Event*) this);
    }
}

VecPlayContinuousX::~VecPlayContinuousX()
{
    delete [] y_;
    delete [] t_;
};

void VecPlayContinuousX::Continuous(double tt) {
    *pd_ = Interpolate(tt);
}

double VecPlayContinuousX::Interpolate(double tt) {
    if (tt >= t_[ubound_index_]) {
        last_index_ = ubound_index_;
        if (last_index_ == 0) {
            return y_[last_index_];
        }
    } else if (tt <= t_[0]) {
        last_index_ = 0;
        return y_[0];
    } else {
        Search(tt);
    }
    double x0 = y_[last_index_ - 1];
    double x1 = y_[last_index_];
    double t0 = t_[last_index_ - 1];
    double t1 = t_[last_index_];
    if (t0 == t1) {
        return (x0 + x1) / 2.;
    }
    return Interp((tt - t0) / (t1 - t0), x0, x1);
}

void VecPlayContinuousX::Search(double tt) {
    while (tt < t_[last_index_]) {
        --last_index_;
    }
    while (tt >= t_[last_index_]) {
        ++last_index_;
    }
}

void VecPlayContinuousX::Deliver(floble_t tt, Branch* branch)
{
    last_index_ = ubound_index_;
    if (discon_indices_) {
        if (discon_index_ < size_) {
            ubound_index_ = (int)discon_indices_[discon_index_++];
            branch->AddEventToQueue(t_[ubound_index_], (Event*) this);
        } else {
            ubound_index_ = size_ - 1;
        }
    } else {
        if (ubound_index_ < size_ - 1) {
            ubound_index_++;
            branch->AddEventToQueue(t_[ubound_index_], (Event*) this);
        }
    }
    Continuous(tt);
}

