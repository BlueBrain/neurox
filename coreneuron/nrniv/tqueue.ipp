/*
Copyright (c) 2016, Blue Brain Project
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef tqueue_ipp_
#define tqueue_ipp_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrniv/tqueue.h"

#if COLLECT_TQueue_STATISTICS
#define STAT(arg) ++arg;
#else
#define STAT(arg) /**/
#endif

// splay tree + bin queue limited to fixed step method
// for event-sets or priority queues
// this starts from the sptqueue.cpp file and adds a bin queue

/* Derived from David Brower's c translation of pascal code by
Douglas Jones.
*/
/* The original c code is included from this file but note that instead
of struct _spblk, we are really using TQItem
*/

template <container C>
TQueue<C>::TQueue() {
    MUTCONSTRUCT(0)
    nshift_ = 0;
    sptree_ = new SPTREE;
    spinit(sptree_);
    binq_ = new BinQ;
    least_ = 0;

#if COLLECT_TQueue_STATISTICS
    nmove = ninsert = nrem = nleast = nbal = ncmplxrem = 0;
    nfastmove = ncompare = nleastsrch = nfind = nfindsrch = 0;
#endif
}

template <container C>
TQueue<C>::~TQueue() {
    SPBLK *q, *q2;
    /// Clear the binq
    for (q = binq_->first(); q; q = q2) {
        q2 = binq_->next(q);
        remove(q);  /// Potentially dereferences freed pointer this->sptree_
    }
    delete binq_;

    /// Clear the splay tree
    while ((q = spdeq(&sptree_->root)) != NULL) {
        delete q;
    }
    delete sptree_;

    /// Clear the priority queue
    while (pq_que_.size()) {
        delete pq_que_.top().second;
        pq_que_.pop();
    }

    MUTDESTRUCT
}

template <container C>
TQItem* TQueue<C>::enqueue_bin(double td, void* d) {
    MUTLOCK
#if COLLECT_TQueue_STATISTICS
    STAT(ninsert);
    record_stat_event(enq, td);
#endif
    TQItem* i = new TQItem;
    i->data_ = d;
    i->t_ = td;
    binq_->enqueue(td, i);
    MUTUNLOCK
    return i;
}

#if COLLECT_TQueue_STATISTICS
template <container C>
void TQueue<C>::record_stat_event(int type, double time) {
    if (time_map_events[type].find(time) == time_map_events[type].end())
        time_map_events[type][time] = 1;
    else
        ++time_map_events[type][time];
}

template <container C>
void TQueue<C>::statistics() {
    printf("insertions=%lu  moves=%lu removals=%lu calls to least=%lu\n", ninsert, nmove, nrem,
           nleast);
    printf("calls to find=%lu\n", nfind);
    printf("comparisons=%d\n", sptree_->enqcmps);
}
#endif

/// Splay tree priority queue implementation
template <>
inline void TQueue<spltree>::move_least_nolock(double tnew) {
    TQItem* b = least();
    if (b) {
        b->t_ = tnew;
        TQItem* nl;
        nl = sphead(sptree_);
        if (nl && (tnew > nl->t_)) {
            least_ = spdeq(&sptree_->root);
            spenq(b, sptree_);
        }
    }
}

/// STL priority queue implementation
template <>
inline void TQueue<pq_que>::move_least_nolock(double tnew) {
    TQItem* b = least();
    if (b) {
        b->t_ = tnew;
        TQItem* nl;
        nl = pq_que_.top().second;
        if (nl && (tnew > nl->t_)) {
            least_ = nl;
            pq_que_.pop();
            pq_que_.push(make_TQPair(b));
        }
    }
}

/// Splay tree priority queue implementation
template <>
inline void TQueue<spltree>::move(TQItem* i, double tnew) {
    MUTLOCK
    STAT(nmove)
    if (i == least_) {
        move_least_nolock(tnew);
    } else if (tnew < least_->t_) {
        spdelete(i, sptree_);
        i->t_ = tnew;
        spenq(least_, sptree_);
        least_ = i;
    } else {
        spdelete(i, sptree_);
        i->t_ = tnew;
        spenq(i, sptree_);
    }
    MUTUNLOCK
}

/// STL priority queue implementation
template <>
inline void TQueue<pq_que>::move(TQItem* i, double tnew) {
    MUTLOCK
    STAT(nmove)
    if (i == least_) {
        move_least_nolock(tnew);
    } else if (tnew < least_->t_) {
        TQItem* qmove = new TQItem;
        qmove->data_ = i->data_;
        qmove->t_ = tnew;
        qmove->cnt_ = i->cnt_;
        i->t_ = -1.;
        pq_que_.push(make_TQPair(least_));
        least_ = qmove;
    } else {
        TQItem* qmove = new TQItem;
        qmove->data_ = i->data_;
        qmove->t_ = tnew;
        qmove->cnt_ = i->cnt_;
        i->t_ = -1.;
        pq_que_.push(make_TQPair(qmove));
    }
    MUTUNLOCK
}

/// Splay tree priority queue implementation
template <>
inline TQItem* TQueue<spltree>::insert(double tt, void* d) {
    MUTLOCK
#if COLLECT_TQueue_STATISTICS
    STAT(ninsert);
    record_stat_event(enq, tt);
#endif
    TQItem* i = new TQItem;
    i->data_ = d;
    i->t_ = tt;
    i->cnt_ = -1;
    if (tt < least_t_nolock()) {
        if (least_) {
            /// Probably storing both time and event which has the time is redundant, but the event
            /// is then returned
            /// to the upper level call stack function. If we were to eliminate i->t_ and i->cnt_
            /// fields,
            /// we need to make sure we are not braking anything.
            spenq(least_, sptree_);
        }
        least_ = i;
    } else {
        spenq(i, sptree_);
    }
    MUTUNLOCK
    return i;
}

/// STL priority queue implementation
template <>
inline TQItem* TQueue<pq_que>::insert(double tt, void* d) {
    MUTLOCK
#if COLLECT_TQueue_STATISTICS
    STAT(ninsert);
    record_stat_event(enq, tt);
#endif
    TQItem* i = new TQItem;
    i->data_ = d;
    i->t_ = tt;
    i->cnt_ = -1;
    if (tt < least_t_nolock()) {
        if (least_) {
            /// Probably storing both time and event which has the time is redundant, but the event
            /// is then returned
            /// to the upper level call stack function. If we were to eliminate i->t_ and i->cnt_
            /// fields,
            /// we need to make sure we are not braking anything.
            pq_que_.push(make_TQPair(least_));
        }
        least_ = i;
    } else {
        pq_que_.push(make_TQPair(i));
    }
    MUTUNLOCK
    return i;
}

/// Splay tree priority queue implementation
template <>
inline void TQueue<spltree>::remove(TQItem* q) {
    MUTLOCK
#if COLLECT_TQueue_STATISTICS
    STAT(nrem);
    record_stat_event(deq, q->t_);
#endif
    if (q) {
        if (q == least_) {
            if (sptree_->root) {
                least_ = spdeq(&sptree_->root);
            } else {
                least_ = NULL;
            }
        } else {
            spdelete(q, sptree_);
        }
        delete q;
    }
    MUTUNLOCK
}

/// STL priority queue implementation
template <>
inline void TQueue<pq_que>::remove(TQItem* q) {
    MUTLOCK
#if COLLECT_TQueue_STATISTICS
    STAT(nrem);
    record_stat_event(deq, q->t_);
#endif
    if (q) {
        if (q == least_) {
            if (pq_que_.size()) {
                least_ = pq_que_.top().second;
                pq_que_.pop();
            } else {
                least_ = NULL;
            }
        } else {
            q->t_ = -1.;
        }
    }
    MUTUNLOCK
}

/// Splay tree priority queue implementation
template <>
inline TQItem* TQueue<spltree>::atomic_dq(double tt) {
    TQItem* q = 0;
    MUTLOCK
    if (least_ && least_->t_ <= tt) {
        q = least_;
#if COLLECT_TQueue_STATISTICS
        STAT(nrem);
        record_stat_event(deq, tt);
#endif
        if (sptree_->root) {
            least_ = spdeq(&sptree_->root);
        } else {
            least_ = NULL;
        }
    }
    MUTUNLOCK
    return q;
}

/// STL priority queue implementation
template <>
inline TQItem* TQueue<pq_que>::atomic_dq(double tt) {
    TQItem* q = 0;
    MUTLOCK
    if (least_ && least_->t_ <= tt) {
        q = least_;
#if COLLECT_TQueue_STATISTICS
        STAT(nrem);
        record_stat_event(deq, tt);
#endif
        //        int qsize = pq_que_.size();
        //        printf("map size: %d\n", msize);
        /// This while loop is to delete events whose times have been moved with the ::move
        /// function,
        /// but in fact events were left in the queue since the only function available is pop
        while (pq_que_.size() && pq_que_.top().second->t_ < 0.) {
            delete pq_que_.top().second;
            pq_que_.pop();
        }
        if (pq_que_.size()) {
            least_ = pq_que_.top().second;
            pq_que_.pop();
        } else {
            least_ = NULL;
        }
    }
    MUTUNLOCK
    return q;
}

#endif
