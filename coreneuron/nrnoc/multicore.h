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

#ifndef multicore_h
#define multicore_h

#include "coreneuron/nrnconf.h"
#include "coreneuron/nrnoc/membfunc.h"
#include "coreneuron/nrnoc/nrnoc_nt.h"

#if defined(__cplusplus)
class NetCon;
class PreSyn;
extern "C" {
#else
typedef void NetCon;
typedef void PreSyn;
#endif

/*
   Point_process._presyn, used only if its NET_RECEIVE sends a net_event, is
   eliminated. Needed only by net_event function. Replaced by
   PreSyn* = nt->presyns + nt->pnt2presyn_ix[pnttype2presyn[pnt->_type]][pnt->_i_instance];
*/
extern int nrn_has_net_event_cnt_; /* how many net_event sender types are there? */
extern int* nrn_has_net_event_; /* the types that send a net_event */
extern int* pnttype2presyn; /* from the type, which array of pnt2presyn_ix are we talking about. */

extern void nrn_threads_create(int n, int parallel);
extern int nrn_nthread;
extern NrnThread* nrn_threads;
extern void nrn_multithread_job(void*(*)(NrnThread*));
extern void nrn_thread_table_check(void);

extern void nrn_threads_free(void);

extern int n_memb_func;
extern Memb_func* memb_func;

#if defined(__cplusplus)
}
#endif

#endif
