/*
Copyright (c) 2014 EPFL-BBP, All rights reserved.

THIS SOFTWARE IS PROVIDED BY THE BLUE BRAIN PROJECT "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE BLUE BRAIN PROJECT
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef netcvode_h
#define netcvode_h

#include "coreneuron/nrniv/sptbinq.h"

#define PRINT_EVENT 1

class DiscreteEvent;
class SelfEventPool;
class NetCvode;
struct InterThreadEvent;

class NetCvodeThreadData {
public:
    int unreffed_event_cnt_;
    TQueue* tqe_;
    std::vector<InterThreadEvent*> inter_thread_events_;
    MUTDEC

    NetCvodeThreadData();
    virtual ~NetCvodeThreadData();
    void interthread_send(double, DiscreteEvent*, NrnThread*);
    void enqueue(NetCvode*, NrnThread*);
};

class NetCvode {
public:
    int print_event_;
    int pcnt_;
    int enqueueing_;
    NetCvodeThreadData* p;
    static double eps_;

    NetCvode(void);
    virtual ~NetCvode();
    void p_construct(int);
    void check_thresh(NrnThread*);
    static double eps(double x) { return eps_*fabs(x); }
    TQItem* event(double tdeliver, DiscreteEvent*, NrnThread*);
    void move_event(TQItem*, double, NrnThread*);
    TQItem* bin_event(double tdeliver, DiscreteEvent*, NrnThread*);
    void deliver_net_events(NrnThread*); // for default staggered time step method
    void deliver_events(double til, NrnThread*); // for initialization events
    bool deliver_event(double til, NrnThread*); //uses TQueue atomically
    void clear_events();
    void init_events();
    void point_receive(int, Point_process*, double*, double);
};

#endif
