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

#include "coreneuron/coreneuron.h"

#if defined(_OPENACC)
#define _PRAGMA_FOR_INIT_ACC_LOOP_ \
    _Pragma("acc parallel loop present(vdata[0:_cntml_padded*nparm]) if(_nt->compute_gpu)")
#define _PRAGMA_FOR_CUR_ACC_LOOP_ \
    _Pragma(                      \
        "acc parallel loop present(vdata[0:_cntml_padded*nparm], ni[0:_cntml_actual], _vec_rhs[0:_nt->end]) if(_nt->compute_gpu) async(stream_id)")
#define _PRAGMA_FOR_JACOB_ACC_LOOP_ \
    _Pragma(                        \
        "acc parallel loop present(vdata[0:_cntml_padded*nparm], ni[0:_cntml_actual], _vec_d[0:_nt->end]) if(_nt->compute_gpu) async(stream_id)")
#else
#define _PRAGMA_FOR_INIT_ACC_LOOP_ _Pragma("")
#define _PRAGMA_FOR_CUR_ACC_LOOP_ _Pragma("")
#define _PRAGMA_FOR_JACOB_ACC_LOOP_ _Pragma("")
#endif

#if !defined(LAYOUT)
/* 1 means AoS, >1 means AoSoA, <= 0 means SOA */
#define LAYOUT 1
#endif
#if LAYOUT >= 1
#define _STRIDE LAYOUT
#else
#define _STRIDE _cntml_padded + _iml
#endif

static const char* mechanism[] = {"0", "capacitance", "cm", 0, "i_cap", 0, 0};
void nrn_alloc_capacitance(double*, Datum*, int);
void nrn_init_capacitance(NrnThread*, Memb_list*, int);
void nrn_jacob_capacitance(NrnThread*, Memb_list*, int);

#define nparm 2

void capacitance_reg(void) {
    int mechtype;
    /* all methods deal with capacitance in special ways */
    register_mech(mechanism, nrn_alloc_capacitance, (mod_f_t)0, (mod_f_t)0, (mod_f_t)0,
                  (mod_f_t)nrn_init_capacitance, -1, 1);
    mechtype = nrn_get_mechtype(mechanism[1]);
    _nrn_layout_reg(mechtype, LAYOUT);
    hoc_register_prop_size(mechtype, nparm, 0);
}

#define cm vdata[0 * _STRIDE]
#define i_cap vdata[1 * _STRIDE]

/*
cj is analogous to 1/dt for cvode and daspk
for fixed step second order it is 2/dt and
for pure implicit fixed step it is 1/dt
It used to be static but is now a thread data variable
*/

void nrn_jacob_capacitance(NrnThread* _nt, Memb_list* ml, int type) {
    (void)type;
    int _cntml_actual = ml->nodecount;
    int _cntml_padded = ml->_nodecount_padded;
    int _iml;
    double* vdata;
    double cfac = .001 * _nt->cj;
    (void)_cntml_padded; /* unused when layout=1*/

    double* _vec_d = _nt->_actual_d;
#if defined(_OPENACC)
    int stream_id = _nt->stream_id;
#endif

    { /*if (use_cachevec) {*/
        int* ni = ml->nodeindices;

#if LAYOUT == 1 /*AoS*/
        for (_iml = 0; _iml < _cntml_actual; _iml++) {
            vdata = ml->data + _iml * nparm;
#else
        vdata = ml->data;
        _PRAGMA_FOR_JACOB_ACC_LOOP_
        for (_iml = 0; _iml < _cntml_actual; _iml++) {
#endif
            _vec_d[ni[_iml]] += cfac * cm;
        }
    }
}

void nrn_init_capacitance(NrnThread* _nt, Memb_list* ml, int type) {
    (void)type;
    int _cntml_actual = ml->nodecount;
    int _cntml_padded = ml->_nodecount_padded;
    int _iml;
    double* vdata;
    (void)_nt;
    (void)type;
    (void)_cntml_padded; /* unused */

#if LAYOUT == 1 /*AoS*/
    for (_iml = 0; _iml < _cntml_actual; _iml++) {
        vdata = ml->data + _iml * nparm;
#else
    vdata = ml->data;
    _PRAGMA_FOR_INIT_ACC_LOOP_
    for (_iml = 0; _iml < _cntml_actual; _iml++) {
#endif
        i_cap = 0;
    }
}

void nrn_cur_capacitance(NrnThread* _nt, Memb_list* ml, int type) {
    (void)type;
    int _cntml_actual = ml->nodecount;
    int _cntml_padded = ml->_nodecount_padded;
    int _iml;
    double* vdata;
    double cfac = .001 * _nt->cj;

    /*@todo: verify cfac is being copied !! */

    (void)_cntml_padded; /* unused when layout=1*/

    /* since rhs is dvm for a full or half implicit step */
    /* (nrn_update_2d() replaces dvi by dvi-dvx) */
    /* no need to distinguish secondorder */
    int* ni = ml->nodeindices;
    double* _vec_rhs = _nt->_actual_rhs;
#if defined(_OPENACC)
    int stream_id = _nt->stream_id;
#endif

#if LAYOUT == 1 /*AoS*/
    for (_iml = 0; _iml < _cntml_actual; _iml++) {
        vdata = ml->data + _iml * nparm;
#else
    vdata = ml->data;
    _PRAGMA_FOR_CUR_ACC_LOOP_
    for (_iml = 0; _iml < _cntml_actual; _iml++) {
#endif
        i_cap = cfac * cm * _vec_rhs[ni[_iml]];
    }
}

/* the rest can be constructed automatically from the above info*/

void nrn_alloc_capacitance(double* data, Datum* pdata, int type) {
    (void)pdata;
    (void)type;       /* unused */
    data[0] = DEF_cm; /*default capacitance/cm^2*/
}
