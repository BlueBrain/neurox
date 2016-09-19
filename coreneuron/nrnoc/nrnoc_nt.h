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

#ifndef nrnoc_nt_h
#define nrnoc_nt_h

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include <stdint.h>

#if defined(__cplusplus)
class NetCon;
class PreSyn;
extern "C" {
#else
typedef void NetCon;
typedef void PreSyn;
#endif

typedef struct NrnThreadMembList{ /* patterned after CvMembList in cvodeobj.h */
    struct NrnThreadMembList* next;
    struct Memb_list* ml;
    int index;
    int *dependencies; /* list of mechanism types that this mechanism depends on*/
    int ndependencies; /* for scheduling we need to know the dependency count */
} NrnThreadMembList;

typedef struct NrnThreadBAList {
    struct Memb_list* ml; /* an item in the NrnThreadMembList */
    struct BAMech* bam;
    struct NrnThreadBAList* next;
} NrnThreadBAList;

/* for OpenACC, in order to avoid an error while update PreSyn, with virtual base
 * class, we are adding helper with flag variable which could be updated on GPU
 */
typedef struct PreSynHelper {
    int flag_;
} PreSynHelper;

/**
 * The NrnThread class
 * Describes all metadata for a computation job, composed by a group of neurons;
 */
typedef struct NrnThread {
    double _t;            ///> current execution time
    double _dt;           ///> user-defined fixed time step
    double cj;            ///> 2.0/dt for secondorder executions, otherwise 1.0/dt

    NrnThreadMembList* tml;
    Memb_list** _ml_list;
    Point_process* pntprocs; // synapses and artificial cells with and without gid
    PreSyn* presyns; // all the output PreSyn with and without gid
    PreSynHelper* presyns_helper;
    int** pnt2presyn_ix; // eliminates Point_process._presyn used only by net_event sender.
    NetCon* netcons;
    double* weights;      ///> array of arguments for netcons;

    int n_pntproc, n_presyn, n_input_presyn, n_netcon, n_weight; // only for model_size

    int ncell;            ///> number of cells (described by the first compartments)
    int end;              ///> number of compartments
    int id;               ///> if of current NrnThread in nrn_threads array;
    int _stop_stepping;

    size_t _ndata;        ///> size of _data array;
    size_t _nidata;       ///> size of _idata array;
    size_t _nvdata;       ///> size of _vdata array;
    size_t n_vecplay;     ///> size of _vecplay array;

    double* _data;        ///> all double data for the branch (RHS, D, A, B, V, Area, and mechanisms)
    int* _idata;          ///> all the Datum to ints index into here
    void** _vdata;        ///> all the Datum to pointers index into here
    void** _vecplay;      ///> array of instances of VecPlayContinuous

    double* _actual_rhs;  ///> right-hand side (solution vector) of Linear Algebra solver
    double* _actual_d;    ///> main diagonal of Linear Algebra spart tridiagonal matrix
    double* _actual_a;    ///> top diagonal of Linear Algebra sparse tridiagonal matrix
    double* _actual_b;    ///> bottom diagonal of Linear Algebra sparse tridiagonal matrix
    double* _actual_v;    ///> diagonal of Linear Algebra sparse tridiagonal matrix
    double* _actual_area; ///> current area per compartment
    double* _shadow_rhs;  ///> Not pointer into _data. Avoid race for multiple POINT_PROCESS in same compartment
    double* _shadow_d;    ///> Not pointer into _data. Avoid race for multiple POINT_PROCESS in same compartment
    int* _v_parent_index; ///> index of parents compartments (if multiSpliX is 0) or NULL (if multiSpliX is 1)
    int* _permute;
    char* _sp13mat;       ///> handle to general sparse matrix
    struct Memb_list* _ecell_memb_list; /* normally nil */

    double _ctime; /* computation time in seconds (using nrnmpi_wtime) */

    NrnThreadBAList* tbl[BEFORE_AFTER_SIZE]; ///> Array of before-after function pointers

    int shadow_rhs_cnt;   ///> added to facilitate the NrnThread transfer to GPU
    int compute_gpu;      ///> define whether to compute with gpus
    int stream_id;        ///> define where the kernel will be launched on GPU stream
    int _net_send_buffer_size;
    int _net_send_buffer_cnt;
    int* _net_send_buffer;
    void* mapping;        ///> section to segment mapping information

} NrnThread;

#if defined(__cplusplus)
}
#endif

#endif
