#include "neurox/neurox.h"
#include <algorithm>
#include <cstring>

using namespace neurox;

VecPlayContinuousX::VecPlayContinuousX(double* pd, IvocVect* yvec, IvocVect* tvec, IvocVect* discon)
    :pd_(pd), y_(yvec), t_(tvec), discon_indices_(discon), last_index_(0), discon_index_(0), ubound_index_(0)
{}

void VecPlayContinuousX::play_init(Branch * branch) {
    last_index_ = 0;
    discon_index_ = 0;
    if (discon_indices_) {
        if (discon_indices_->size() > 0) {
            ubound_index_ = (int)(*discon_indices_)[discon_index_++];
            branch->addEventToQueue( (*t_)[ubound_index_], (Event*) this);
        } else {
            ubound_index_ = t_->size() - 1;
        }
    } else {
        ubound_index_ = 0;
        branch->addEventToQueue((*t_)[ubound_index_], (Event*) this);
    }
}

VecPlayContinuousX::~VecPlayContinuousX()
{};

void VecPlayContinuousX::continuous(double tt) {
    *pd_ = interpolate(tt);
}

double VecPlayContinuousX::interpolate(double tt) {
    if (tt >= (*t_)[ubound_index_]) {
        last_index_ = ubound_index_;
        if (last_index_ == 0) {
            return (*y_)[last_index_];
        }
    } else if (tt <= (*t_)[0]) {
        last_index_ = 0;
        return (*y_)[0];
    } else {
        search(tt);
    }
    double x0 = (*y_)[last_index_ - 1];
    double x1 = (*y_)[last_index_];
    double t0 = (*t_)[last_index_ - 1];
    double t1 = (*t_)[last_index_];
    if (t0 == t1) {
        return (x0 + x1) / 2.;
    }
    return interp((tt - t0) / (t1 - t0), x0, x1);
}

void VecPlayContinuousX::search(double tt) {
    while (tt < (*t_)[last_index_]) {
        --last_index_;
    }
    while (tt >= (*t_)[last_index_]) {
        ++last_index_;
    }
}

void VecPlayContinuousX::deliver(floble_t tt, Branch* branch)
{
    last_index_ = ubound_index_;
    if (discon_indices_) {
        if (discon_index_ < discon_indices_->size()) {
            ubound_index_ = (int)(*discon_indices_)[discon_index_++];
            branch->addEventToQueue((*t_)[ubound_index_], (Event*) this);
        } else {
            ubound_index_ = t_->size() - 1;
        }
    } else {
        if (ubound_index_ < t_->size() - 1) {
            ubound_index_++;
            branch->addEventToQueue((*t_)[ubound_index_], (Event*) this);
        }
    }
    continuous(tt);
}

