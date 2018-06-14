/* Created by Language version: 6.2.0 */
/* VECTORIZED */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#undef PI
 
#include "coreneuron/utils/randoms/nrnran123.h"
#include "coreneuron/nrnoc/md1redef.h"
#include "coreneuron/nrnconf.h"
#include "coreneuron/nrnoc/membfunc.h"
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrniv/nrniv_decl.h"
#include "coreneuron/nrniv/ivocvect.h"
#include "coreneuron/nrniv/nrn_acc_manager.h"
#include "coreneuron/mech/cfile/scoplib.h"

#include "coreneuron/scopmath_core/newton_struct.h"
#include "coreneuron/nrnoc/md2redef.h"
#include "coreneuron/nrnoc/register_mech.hpp"
#include "_kinderiv.h"
#if !NRNGPU
#if !defined(DISABLE_HOC_EXP)
#undef exp
#define exp hoc_Exp
#endif
#endif
 #define _threadargsproto_namespace int _iml, int _cntml_padded, double* _p, coreneuron::Datum* _ppvar, coreneuron::ThreadDatum* _thread, coreneuron::NrnThread* _nt, double v
 static void bbcore_read(double *, int*, int*, int*, _threadargsproto_namespace);
 static void bbcore_write(double *, int*, int*, int*, _threadargsproto_namespace);
 namespace coreneuron {
 
#define _thread_present_ /**/ , _slist1[0:2], _dlist1[0:2] 
 
#if defined(_OPENACC) && !defined(DISABLE_OPENACC)
#include <openacc.h>
#if defined(PG_ACC_BUGS)
#define _PRAGMA_FOR_INIT_ACC_LOOP_ _Pragma("acc parallel loop present(_ni[0:_cntml_actual], _nt_data[0:_nt->_ndata], _p[0:_cntml_padded*_psize], _ppvar[0:_cntml_padded*_ppsize], _vec_v[0:_nt->end], nrn_ion_global_map[0:nrn_ion_global_map_size][0:3], _nt[0:1] _thread_present_) if(_nt->compute_gpu)")
#else
#define _PRAGMA_FOR_INIT_ACC_LOOP_ _Pragma("acc parallel loop present(_ni[0:_cntml_actual], _nt_data[0:_nt->_ndata], _p[0:_cntml_padded*_psize], _ppvar[0:_cntml_padded*_ppsize], _vec_v[0:_nt->end], nrn_ion_global_map[0:nrn_ion_global_map_size], _nt[0:1] _thread_present_) if(_nt->compute_gpu)")
#endif
#define _PRAGMA_FOR_STATE_ACC_LOOP_ _Pragma("acc parallel loop present(_ni[0:_cntml_actual], _nt_data[0:_nt->_ndata], _p[0:_cntml_padded*_psize], _ppvar[0:_cntml_padded*_ppsize], _vec_v[0:_nt->end], _nt[0:1], _ml[0:1] _thread_present_) if(_nt->compute_gpu) async(stream_id)")
#define _PRAGMA_FOR_CUR_ACC_LOOP_ _Pragma("acc parallel loop present(_ni[0:_cntml_actual], _nt_data[0:_nt->_ndata], _p[0:_cntml_padded*_psize], _ppvar[0:_cntml_padded*_ppsize], _vec_v[0:_nt->end], _vec_d[0:_nt->end], _vec_rhs[0:_nt->end], _nt[0:1] _thread_present_) if(_nt->compute_gpu) async(stream_id)")
#define _PRAGMA_FOR_CUR_SYN_ACC_LOOP_ _Pragma("acc parallel loop present(_ni[0:_cntml_actual], _nt_data[0:_nt->_ndata], _p[0:_cntml_padded*_psize], _ppvar[0:_cntml_padded*_ppsize], _vec_v[0:_nt->end], _vec_shadow_rhs[0:_nt->shadow_rhs_cnt], _vec_shadow_d[0:_nt->shadow_rhs_cnt], _vec_d[0:_nt->end], _vec_rhs[0:_nt->end], _nt[0:1]) if(_nt->compute_gpu) async(stream_id)")
#define _PRAGMA_FOR_NETRECV_ACC_LOOP_ _Pragma("acc parallel loop present(_pnt[0:_pnt_length], _nrb[0:1], _nt[0:1], nrn_threads[0:nrn_nthread]) if(_nt->compute_gpu) async(stream_id)")
#define _ACC_GLOBALS_UPDATE_ if (_nt->compute_gpu) {_acc_globals_update();}
#else
#define _PRAGMA_FOR_INIT_ACC_LOOP_ _Pragma("")
#define _PRAGMA_FOR_STATE_ACC_LOOP_ _Pragma("")
#define _PRAGMA_FOR_CUR_ACC_LOOP_ _Pragma("")
#define _PRAGMA_FOR_CUR_SYN_ACC_LOOP_ _Pragma("")
#define _PRAGMA_FOR_NETRECV_ACC_LOOP_ _Pragma("")
#define _ACC_GLOBALS_UPDATE_ ;
#endif
 
#if defined(__clang__)
#define _PRAGMA_FOR_VECTOR_LOOP_ _Pragma("clang loop vectorize(enable)")
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#define _PRAGMA_FOR_VECTOR_LOOP_ _Pragma("ivdep")
#elif defined(__IBMC__) || defined(__IBMCPP__)
#define _PRAGMA_FOR_VECTOR_LOOP_ _Pragma("ibm independent_loop")
#elif defined(__PGI)
#define _PRAGMA_FOR_VECTOR_LOOP_ _Pragma("vector")
#elif defined(_CRAYC)
#define _PRAGMA_FOR_VECTOR_LOOP_ _Pragma("_CRI ivdep")
#elif defined(__GNUC__) || defined(__GNUG__)
#define _PRAGMA_FOR_VECTOR_LOOP_ _Pragma("GCC ivdep")
#else
#define _PRAGMA_FOR_VECTOR_LOOP_
#endif // _PRAGMA_FOR_VECTOR_LOOP_
 
#if !defined(LAYOUT)
/* 1 means AoS, >1 means AoSoA, <= 0 means SOA */
#define LAYOUT 1
#endif
#if LAYOUT >= 1
#define _STRIDE LAYOUT
#else
#define _STRIDE _cntml_padded + _iml
#endif
 
#define nrn_init _nrn_init__StochKv3
#define nrn_cur _nrn_cur__StochKv3
#define _nrn_current _nrn_current__StochKv3
#define nrn_jacob _nrn_jacob__StochKv3
#define nrn_state _nrn_state__StochKv3
#define initmodel initmodel__StochKv3
#define _net_receive _net_receive__StochKv3
#define nrn_state_launcher nrn_state_StochKv3_launcher
#define nrn_cur_launcher nrn_cur_StochKv3_launcher
#define nrn_jacob_launcher nrn_jacob_StochKv3_launcher 
#define ChkProb ChkProb_StochKv3 
#define _f_trates _f_trates_StochKv3 
#define setRNG setRNG_StochKv3 
#define states states_StochKv3 
#define trates trates_StochKv3 
#define _ode_matsol1 _nrn_ode_matsol1__StochKv3
#define _ode_spec1 _nrn_ode_spec1__StochKv3

#undef _threadargscomma_
#undef _threadargsprotocomma_
#undef _threadargs_
#undef _threadargsproto_
 
#define _threadargscomma_ _iml, _cntml_padded, _p, _ppvar, _thread, _nt, v,
#define _threadargsprotocomma_ int _iml, int _cntml_padded, double* _p, Datum* _ppvar, ThreadDatum* _thread, NrnThread* _nt, double v,
#define _threadargs_ _iml, _cntml_padded, _p, _ppvar, _thread, _nt, v
#define _threadargsproto_ int _iml, int _cntml_padded, double* _p, Datum* _ppvar, ThreadDatum* _thread, NrnThread* _nt, double v
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gamma _p[0*_STRIDE]
#define eta _p[1*_STRIDE]
#define gkbar _p[2*_STRIDE]
#define deterministic _p[3*_STRIDE]
#define an _p[4*_STRIDE]
#define bn _p[5*_STRIDE]
#define al _p[6*_STRIDE]
#define bl _p[7*_STRIDE]
#define ik _p[8*_STRIDE]
#define gk _p[9*_STRIDE]
#define ninf _p[10*_STRIDE]
#define ntau _p[11*_STRIDE]
#define linf _p[12*_STRIDE]
#define ltau _p[13*_STRIDE]
#define N _p[14*_STRIDE]
#define P_an _p[15*_STRIDE]
#define P_bn _p[16*_STRIDE]
#define P_al _p[17*_STRIDE]
#define P_bl _p[18*_STRIDE]
#define n _p[19*_STRIDE]
#define l _p[20*_STRIDE]
#define ek _p[21*_STRIDE]
#define scale_dens _p[22*_STRIDE]
#define usingR123 _p[23*_STRIDE]
#define Dn _p[24*_STRIDE]
#define Dl _p[25*_STRIDE]
#define N0L0 _p[26*_STRIDE]
#define N1L0 _p[27*_STRIDE]
#define N0L1 _p[28*_STRIDE]
#define N1L1 _p[29*_STRIDE]
#define n0l0_n1l0 _p[30*_STRIDE]
#define n0l0_n0l1 _p[31*_STRIDE]
#define n1l0_n1l1 _p[32*_STRIDE]
#define n1l0_n0l0 _p[33*_STRIDE]
#define n0l1_n1l1 _p[34*_STRIDE]
#define n0l1_n0l0 _p[35*_STRIDE]
#define n1l1_n0l1 _p[36*_STRIDE]
#define n1l1_n1l0 _p[37*_STRIDE]
#define _v_unused _p[38*_STRIDE]
#define _g_unused _p[39*_STRIDE]
 
#ifndef NRN_PRCELLSTATE
#define NRN_PRCELLSTATE 0
#endif
#if NRN_PRCELLSTATE
#define _PRCELLSTATE_V _v_unused = _v;
#define _PRCELLSTATE_G _g_unused = _g;
#else
#define _PRCELLSTATE_V /**/
#define _PRCELLSTATE_G /**/
#endif
#define _ion_ek		_nt_data[_ppvar[0*_STRIDE]]
#define _ion_ik	_nt_data[_ppvar[1*_STRIDE]]
#define _ion_dikdv	_nt_data[_ppvar[2*_STRIDE]]
#define _p_rng	_nt->_vdata[_ppvar[3*_STRIDE]]
#define area	_nt->_data[_ppvar[4*_STRIDE]]
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 static int hoc_nrnpointerindex =  3;
 static ThreadDatum* _extcall_thread;
 /* external NEURON variables */
 
#if 0 /*BBCORE*/
 /* declaration of user functions */
 static void _hoc_BnlDev(void);
 static void _hoc_ChkProb(void);
 static void _hoc_bbsavestate(void);
 static void _hoc_brand(void);
 static void _hoc_setRNG(void);
 static void _hoc_strap(void);
 static void _hoc_trates(void);
 static void _hoc_urand(void);
 
#endif /*BBCORE*/
 static int _mechtype;
 
#if 0 /*BBCORE*/
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_StochKv3", _hoc_setdata,
 "BnlDev_StochKv3", _hoc_BnlDev,
 "ChkProb_StochKv3", _hoc_ChkProb,
 "bbsavestate_StochKv3", _hoc_bbsavestate,
 "brand_StochKv3", _hoc_brand,
 "setRNG_StochKv3", _hoc_setRNG,
 "strap_StochKv3", _hoc_strap,
 "trates_StochKv3", _hoc_trates,
 "urand_StochKv3", _hoc_urand,
 0, 0
};
 
#endif /*BBCORE*/
#define BnlDev BnlDev_StochKv3
#define bbsavestate bbsavestate_StochKv3
#define brand brand_StochKv3
#define strap strap_StochKv3
#define urand urand_StochKv3
 inline double BnlDev( _threadargsprotocomma_ double , double );
 inline double bbsavestate( _threadargsproto_ );
 inline double brand( _threadargsprotocomma_ double , double );
 inline double strap( _threadargsprotocomma_ double );
 inline double urand( _threadargsproto_ );
 
static void _check_trates(_threadargsproto_); 
static void _check_table_thread(int _iml, int _cntml_padded, double* _p, Datum* _ppvar, ThreadDatum* _thread, NrnThread* _nt, int v) {
   _check_trates(_threadargs_);
 }
 /* declare global and static user variables */
#define usetable usetable_StochKv3
 double usetable = 1;
 #pragma acc declare copyin (usetable)
#define vmax vmax_StochKv3
 double vmax = 100;
 #pragma acc declare copyin (vmax)
#define vmin vmin_StochKv3
 double vmin = -120;
 #pragma acc declare copyin (vmin)
 
static void _acc_globals_update() {
 #pragma acc update device (usetable) if(nrn_threads->compute_gpu)
 #pragma acc update device (vmax) if(nrn_threads->compute_gpu)
 #pragma acc update device (vmin) if(nrn_threads->compute_gpu)
 }
 
#if 0 /*BBCORE*/
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_StochKv3", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "vmin_StochKv3", "mV",
 "vmax_StochKv3", "mV",
 "gamma_StochKv3", "pS",
 "eta_StochKv3", "1/um2",
 "gkbar_StochKv3", "S/cm2",
 "an_StochKv3", "/ms",
 "bn_StochKv3", "/ms",
 "al_StochKv3", "/ms",
 "bl_StochKv3", "/ms",
 "ik_StochKv3", "mA/cm2",
 "gk_StochKv3", "S/cm2",
 "ntau_StochKv3", "ms",
 "ltau_StochKv3", "ms",
 0,0
};
 
#endif /*BBCORE*/
 static double delta_t = 1;
 static double l0 = 0;
 static double n0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "vmin_StochKv3", &vmin_StochKv3,
 "vmax_StochKv3", &vmax_StochKv3,
 "usetable_StochKv3", &usetable_StochKv3,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(double*, Datum*, int);
void nrn_init(NrnThread*, Memb_list*, int);
void nrn_state(NrnThread*, Memb_list*, int);
 void nrn_cur(NrnThread*, Memb_list*, int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "6.2.0",
"StochKv3",
 "gamma_StochKv3",
 "eta_StochKv3",
 "gkbar_StochKv3",
 "deterministic_StochKv3",
 0,
 "an_StochKv3",
 "bn_StochKv3",
 "al_StochKv3",
 "bl_StochKv3",
 "ik_StochKv3",
 "gk_StochKv3",
 "ninf_StochKv3",
 "ntau_StochKv3",
 "linf_StochKv3",
 "ltau_StochKv3",
 "N_StochKv3",
 "P_an_StochKv3",
 "P_bn_StochKv3",
 "P_al_StochKv3",
 "P_bl_StochKv3",
 0,
 "n_StochKv3",
 "l_StochKv3",
 0,
 "rng_StochKv3",
 0};
 static int _k_type;
 
static void nrn_alloc(double* _p, Datum* _ppvar, int _type) {
 
#if 0 /*BBCORE*/
 	/*initialize range parameters*/
 	gamma = 50;
 	eta = 0;
 	gkbar = 0.01;
 	deterministic = 0;
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
#endif /* BBCORE */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 
#define _psize 40
#define _ppsize 5
 void _StochKv3_reg() {
	int _vectorized = 1;
  _initlists();
 _mechtype = nrn_get_mechtype(_mechanism[1]);
 if (_mechtype == -1) return;
 _nrn_layout_reg(_mechtype, LAYOUT);
 _k_type = nrn_get_mechtype("k_ion"); 
#if 0 /*BBCORE*/
 	ion_reg("k", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 
#endif /*BBCORE*/
 	register_mech(_mechanism, nrn_alloc,nrn_cur, NULL, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
     _nrn_thread_table_reg(_mechtype, _check_table_thread);
   hoc_reg_bbcore_read(_mechtype, bbcore_read);
   hoc_reg_bbcore_write(_mechtype, bbcore_write);
  hoc_register_prop_size(_mechtype, _psize, _ppsize);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "bbcorepointer");
  hoc_register_dparam_semantics(_mechtype, 4, "area");
 	hoc_register_var(hoc_scdoub, hoc_vdoub, NULL);
 }
 static double *_t_ntau;
 static double *_t_ltau;
 static double *_t_ninf;
 static double *_t_linf;
 static double *_t_al;
 static double *_t_bl;
 static double *_t_an;
 static double *_t_bn;
static char *modelname = "StochKv3.mod";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static inline int ChkProb(_threadargsprotocomma_ double);
static inline int _f_trates(_threadargsprotocomma_ double);
static inline int setRNG(_threadargsproto_);
static inline int trates(_threadargsprotocomma_ double);
 
int _ode_spec1(_threadargsproto_);
/*int _ode_matsol1(_threadargsproto_);*/
 static void _n_trates(_threadargsprotocomma_ double _lv);
 
#define _slist1 _slist1_StochKv3
int* _slist1;
#pragma acc declare create(_slist1)

#define _dlist1 _dlist1_StochKv3
int* _dlist1;
#pragma acc declare create(_dlist1)
 static inline int states(_threadargsproto_);
 } using namespace coreneuron; 
/*VERBATIM*/
#include "nrnran123.h"
extern int cvode_active_;

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef CORENEURON_BUILD
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
#endif

 namespace coreneuron { 
/*CVODE*/
 int _ode_spec1 (_threadargsproto_) {int _reset = 0; {
   trates ( _threadargscomma_ v ) ;
   Dl = al - ( al + bl ) * l ;
   Dn = an - ( an + bn ) * n ;
   if ( deterministic  || dt > 1.0 ) {
     N1L1 = n * l * N ;
     N1L0 = n * ( 1.0 - l ) * N ;
     N0L1 = ( 1.0 - n ) * l * N ;
     }
   else {
     N1L1 = floor ( N1L1 + 0.5 ) ;
     N1L0 = floor ( N1L0 + 0.5 ) ;
     N0L1 = floor ( N0L1 + 0.5 ) ;
     N0L0 = N - N1L1 - N1L0 - N0L1 ;
     P_an = strap ( _threadargscomma_ an * dt ) ;
     P_bn = strap ( _threadargscomma_ bn * dt ) ;
     ChkProb ( _threadargscomma_ P_an ) ;
     ChkProb ( _threadargscomma_ P_bn ) ;
     n0l0_n1l0 = BnlDev ( _threadargscomma_ P_an , N0L0 ) ;
     n0l1_n1l1 = BnlDev ( _threadargscomma_ P_an , N0L1 ) ;
     n1l1_n0l1 = BnlDev ( _threadargscomma_ P_bn , N1L1 ) ;
     n1l0_n0l0 = BnlDev ( _threadargscomma_ P_bn , N1L0 ) ;
     N0L0 = strap ( _threadargscomma_ N0L0 - n0l0_n1l0 + n1l0_n0l0 ) ;
     N1L0 = strap ( _threadargscomma_ N1L0 - n1l0_n0l0 + n0l0_n1l0 ) ;
     N0L1 = strap ( _threadargscomma_ N0L1 - n0l1_n1l1 + n1l1_n0l1 ) ;
     N1L1 = strap ( _threadargscomma_ N1L1 - n1l1_n0l1 + n0l1_n1l1 ) ;
     P_al = strap ( _threadargscomma_ al * dt ) ;
     P_bl = strap ( _threadargscomma_ bl * dt ) ;
     ChkProb ( _threadargscomma_ P_al ) ;
     ChkProb ( _threadargscomma_ P_bl ) ;
     n0l0_n0l1 = BnlDev ( _threadargscomma_ P_al , N0L0 - n0l0_n1l0 ) ;
     n1l0_n1l1 = BnlDev ( _threadargscomma_ P_al , N1L0 - n1l0_n0l0 ) ;
     n0l1_n0l0 = BnlDev ( _threadargscomma_ P_bl , N0L1 - n0l1_n1l1 ) ;
     n1l1_n1l0 = BnlDev ( _threadargscomma_ P_bl , N1L1 - n1l1_n0l1 ) ;
     N0L0 = strap ( _threadargscomma_ N0L0 - n0l0_n0l1 + n0l1_n0l0 ) ;
     N1L0 = strap ( _threadargscomma_ N1L0 - n1l0_n1l1 + n1l1_n1l0 ) ;
     N0L1 = strap ( _threadargscomma_ N0L1 - n0l1_n0l0 + n0l0_n0l1 ) ;
     N1L1 = strap ( _threadargscomma_ N1L1 - n1l1_n1l0 + n1l0_n1l1 ) ;
     }
   N0L0 = N - N1L1 - N1L0 - N0L1 ;
   }
 return _reset;
}
 int _ode_matsol1 (_threadargsproto_) {
 trates ( _threadargscomma_ v ) ;
 Dl = Dl  / (1. - dt*( ( - (( al + bl ))*(1.0) ) )) ;
 Dn = Dn  / (1. - dt*( ( - (( an + bn ))*(1.0) ) )) ;
 return 0;
}
 /*END CVODE*/
 static int states (_threadargsproto_) { {
   trates ( _threadargscomma_ v ) ;
    l = l + (1. - exp(dt*(( - (( al + bl ))*(1.0) ))))*(- ( al ) / ( ( - (( al + bl ))*(1.0)) ) - l) ;
    n = n + (1. - exp(dt*(( - (( an + bn ))*(1.0) ))))*(- ( an ) / ( ( - (( an + bn ))*(1.0)) ) - n) ;
   if ( deterministic  || dt > 1.0 ) {
     N1L1 = n * l * N ;
     N1L0 = n * ( 1.0 - l ) * N ;
     N0L1 = ( 1.0 - n ) * l * N ;
     }
   else {
     N1L1 = floor ( N1L1 + 0.5 ) ;
     N1L0 = floor ( N1L0 + 0.5 ) ;
     N0L1 = floor ( N0L1 + 0.5 ) ;
     N0L0 = N - N1L1 - N1L0 - N0L1 ;
     P_an = strap ( _threadargscomma_ an * dt ) ;
     P_bn = strap ( _threadargscomma_ bn * dt ) ;
     ChkProb ( _threadargscomma_ P_an ) ;
     ChkProb ( _threadargscomma_ P_bn ) ;
     n0l0_n1l0 = BnlDev ( _threadargscomma_ P_an , N0L0 ) ;
     n0l1_n1l1 = BnlDev ( _threadargscomma_ P_an , N0L1 ) ;
     n1l1_n0l1 = BnlDev ( _threadargscomma_ P_bn , N1L1 ) ;
     n1l0_n0l0 = BnlDev ( _threadargscomma_ P_bn , N1L0 ) ;
     N0L0 = strap ( _threadargscomma_ N0L0 - n0l0_n1l0 + n1l0_n0l0 ) ;
     N1L0 = strap ( _threadargscomma_ N1L0 - n1l0_n0l0 + n0l0_n1l0 ) ;
     N0L1 = strap ( _threadargscomma_ N0L1 - n0l1_n1l1 + n1l1_n0l1 ) ;
     N1L1 = strap ( _threadargscomma_ N1L1 - n1l1_n0l1 + n0l1_n1l1 ) ;
     P_al = strap ( _threadargscomma_ al * dt ) ;
     P_bl = strap ( _threadargscomma_ bl * dt ) ;
     ChkProb ( _threadargscomma_ P_al ) ;
     ChkProb ( _threadargscomma_ P_bl ) ;
     n0l0_n0l1 = BnlDev ( _threadargscomma_ P_al , N0L0 - n0l0_n1l0 ) ;
     n1l0_n1l1 = BnlDev ( _threadargscomma_ P_al , N1L0 - n1l0_n0l0 ) ;
     n0l1_n0l0 = BnlDev ( _threadargscomma_ P_bl , N0L1 - n0l1_n1l1 ) ;
     n1l1_n1l0 = BnlDev ( _threadargscomma_ P_bl , N1L1 - n1l1_n0l1 ) ;
     N0L0 = strap ( _threadargscomma_ N0L0 - n0l0_n0l1 + n0l1_n0l0 ) ;
     N1L0 = strap ( _threadargscomma_ N1L0 - n1l0_n1l1 + n1l1_n1l0 ) ;
     N0L1 = strap ( _threadargscomma_ N0L1 - n0l1_n0l0 + n0l0_n0l1 ) ;
     N1L1 = strap ( _threadargscomma_ N1L1 - n1l1_n1l0 + n1l0_n1l1 ) ;
     }
   N0L0 = N - N1L1 - N1L0 - N0L1 ;
   }
  return 0;
}
 static double _mfac_trates, _tmin_trates;
  static void _check_trates(_threadargsproto_) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_dt;
  if (!usetable) {return;}
  if (_sav_dt != dt) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_trates =  vmin ;
   _tmax =  vmax ;
   _dx = (_tmax - _tmin_trates)/199.; _mfac_trates = 1./_dx;
   for (_i=0, _x=_tmin_trates; _i < 200; _x += _dx, _i++) {
    _f_trates(_threadargs_, _x);
    _t_ntau[_i] = ntau;
    _t_ltau[_i] = ltau;
    _t_ninf[_i] = ninf;
    _t_linf[_i] = linf;
    _t_al[_i] = al;
    _t_bl[_i] = bl;
    _t_an[_i] = an;
    _t_bn[_i] = bn;
   }
   _sav_dt = dt;
  }
 }

 static int trates(_threadargsproto_, double _lv) { 
#if 0
_check_trates(_threadargs_);
#endif
 _n_trates(_threadargs_, _lv);
 return 0;
 }

 static void _n_trates(_threadargsproto_, double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_trates(_threadargs_, _lv); return; 
}
 _xi = _mfac_trates * (_lv - _tmin_trates);
 if (isnan(_xi)) {
  ntau = _xi;
  ltau = _xi;
  ninf = _xi;
  linf = _xi;
  al = _xi;
  bl = _xi;
  an = _xi;
  bn = _xi;
  return;
 }
 if (_xi <= 0.) {
 ntau = _t_ntau[0];
 ltau = _t_ltau[0];
 ninf = _t_ninf[0];
 linf = _t_linf[0];
 al = _t_al[0];
 bl = _t_bl[0];
 an = _t_an[0];
 bn = _t_bn[0];
 return; }
 if (_xi >= 199.) {
 ntau = _t_ntau[199];
 ltau = _t_ltau[199];
 ninf = _t_ninf[199];
 linf = _t_linf[199];
 al = _t_al[199];
 bl = _t_bl[199];
 an = _t_an[199];
 bn = _t_bn[199];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 ntau = _t_ntau[_i] + _theta*(_t_ntau[_i+1] - _t_ntau[_i]);
 ltau = _t_ltau[_i] + _theta*(_t_ltau[_i+1] - _t_ltau[_i]);
 ninf = _t_ninf[_i] + _theta*(_t_ninf[_i+1] - _t_ninf[_i]);
 linf = _t_linf[_i] + _theta*(_t_linf[_i+1] - _t_linf[_i]);
 al = _t_al[_i] + _theta*(_t_al[_i+1] - _t_al[_i]);
 bl = _t_bl[_i] + _theta*(_t_bl[_i+1] - _t_bl[_i]);
 an = _t_an[_i] + _theta*(_t_an[_i+1] - _t_an[_i]);
 bn = _t_bn[_i] + _theta*(_t_bn[_i+1] - _t_bn[_i]);
 }

 
static int  _f_trates ( _threadargsprotocomma_ double _lv ) {
   _lv = _lv + 10.0 ;
   linf = 1.0 / ( 1.0 + exp ( ( - 30.0 - _lv ) / 10.0 ) ) ;
   ltau = 0.346 * exp ( - _lv / ( 18.272 ) ) + 2.09 ;
   ninf = 1.0 / ( 1.0 + exp ( 0.0878 * ( _lv + 55.1 ) ) ) ;
   ntau = 2.1 * exp ( - _lv / 21.2 ) + 4.627 ;
   _lv = _lv - 10.0 ;
   al = linf / ltau ;
   bl = 1.0 / ltau - al ;
   an = ninf / ntau ;
   bn = 1.0 / ntau - an ;
    return 0; }
 
#if 0 /*BBCORE*/
 
static void _hoc_trates(void) {
  double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 
#if 1
 _check_trates(_threadargs_);
#endif
 _r = 1.;
 trates ( _threadargs_, *getarg(1) ;
 hoc_retpushx(_r);
}
 
#endif /*BBCORE*/
 
double strap ( _threadargsprotocomma_ double _lx ) {
   double _lstrap;
 if ( _lx < 0.0 ) {
     _lstrap = 0.0 ;
     
/*VERBATIM*/
        fprintf (stderr,"skv.mod:strap: negative state");
 }
   else {
     _lstrap = _lx ;
     }
   
return _lstrap;
 }
 
#if 0 /*BBCORE*/
 
static void _hoc_strap(void) {
  double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  strap ( _threadargs_, *getarg(1) ;
 hoc_retpushx(_r);
}
 
#endif /*BBCORE*/
 
static int  ChkProb ( _threadargsprotocomma_ double _lp ) {
   if ( _lp < 0.0  || _lp > 1.0 ) {
     
/*VERBATIM*/
    fprintf(stderr, "StochKv2.mod:ChkProb: argument not a probability.\n");
 }
    return 0; }
 
#if 0 /*BBCORE*/
 
static void _hoc_ChkProb(void) {
  double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 ChkProb ( _threadargs_, *getarg(1) ;
 hoc_retpushx(_r);
}
 
#endif /*BBCORE*/
 
static int  setRNG ( _threadargsproto_ ) {
   
/*VERBATIM*/
    // For compatibility, allow for either MCellRan4 or Random123.  Distinguish by the arg types
    // Object => MCellRan4, seeds (double) => Random123
#ifndef CORENEURON_BUILD
    usingR123 = 0;
    if( ifarg(1) && hoc_is_double_arg(1) ) {
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
        uint32_t a2 = 0;
        uint32_t a3 = 0;

        if (*pv) {
            nrnran123_deletestream(*pv);
            *pv = (nrnran123_State*)0;
        }
        if (ifarg(2)) {
            a2 = (uint32_t)*getarg(2);
        }
        if (ifarg(3)) {
            a3 = (uint32_t)*getarg(3);
        }
        *pv = nrnran123_newstream3((uint32_t)*getarg(1), a2, a3);
        usingR123 = 1;
    } else if( ifarg(1) ) {
        void** pv = (void**)(&_p_rng);
        *pv = nrn_random_arg(1);
    } else {
        void** pv = (void**)(&_p_rng);
        *pv = (void*)0;
    }
#endif
  return 0; }
 
#if 0 /*BBCORE*/
 
static void _hoc_setRNG(void) {
  double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 setRNG ( _threadargs_ ;
 hoc_retpushx(_r);
}
 
#endif /*BBCORE*/
 
double urand ( _threadargsproto_ ) {
   double _lurand;
 
/*VERBATIM*/
    double value;
    if( usingR123 ) {
        value = nrnran123_dblpick((nrnran123_State*)_p_rng);
    } else if (_p_rng) {
#ifndef CORENEURON_BUILD
        value = nrn_random_pick(_p_rng);
#endif
    } else {
        value = 0.5;
    }
    _lurand = value;
 
return _lurand;
 }
 
#if 0 /*BBCORE*/
 
static void _hoc_urand(void) {
  double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  urand ( _threadargs_ ;
 hoc_retpushx(_r);
}
 
#endif /*BBCORE*/
 } using namespace coreneuron; 
/*VERBATIM*/
static void bbcore_write(double* x, int* d, int* xx, int* offset, _threadargsproto_) {
    if (d) {
        uint32_t* di = ((uint32_t*)d) + *offset;
      // temporary just enough to see how much space is being used
      if (!_p_rng) {
        di[0] = 0; di[1] = 0, di[2] = 0;
      }else{
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
        nrnran123_getids3(*pv, di, di+1, di+2);
        // write stream sequence
        unsigned char which;
        nrnran123_getseq(*pv, di+3, &which);
        di[4] = (int)which;
      }
      //printf("StochKv3.mod %p: bbcore_write offset=%d %d %d\n", _p, *offset, d?di[0]:-1, d?di[1]:-1);
    }
    *offset += 5;
}
static void bbcore_read(double* x, int* d, int* xx, int* offset, _threadargsproto_) {
    assert(!_p_rng);
    uint32_t* di = ((uint32_t*)d) + *offset;
        if (di[0] != 0 || di[1] != 0|| di[2] != 0)
        {
      nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
      *pv = nrnran123_newstream3(di[0], di[1], di[2]);
      nrnran123_setseq(*pv, di[3], (char)di[4]);
        }
      //printf("StochKv3.mod %p: bbcore_read offset=%d %d %d\n", _p, *offset, di[0], di[1]);
    *offset += 5;
}
 namespace coreneuron { 
double brand ( _threadargsprotocomma_ double _lP , double _lN ) {
   double _lbrand;
 
/*VERBATIM*/
        /*
        :Supports separate independent but reproducible streams for
        : each instance. However, the corresponding hoc Random
        : distribution MUST be set to Random.uniform(0,1)
        */

        // Should probably be optimized
        double value = 0.0;
        int i;
        for (i = 0; i < _lN; i++) {
           if (urand(_threadargs_) < _lP) {
              value = value + 1;
           }
        }
        return(value);

 _lbrand = value ;
   
return _lbrand;
 }
 
#if 0 /*BBCORE*/
 
static void _hoc_brand(void) {
  double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  brand ( _threadargs_, *getarg(1) , *getarg(2) ;
 hoc_retpushx(_r);
}
 
#endif /*BBCORE*/
 } using namespace coreneuron; 
/*VERBATIM*/
#define        PI 3.141592654
#define        r_ia     16807
#define        r_im     2147483647
#define        r_am     (1.0/r_im)
#define        r_iq     127773
#define        r_ir     2836
#define        r_ntab   32
#define        r_ndiv   (1+(r_im-1)/r_ntab)
#define        r_eps    1.2e-7
#define        r_rnmx   (1.0-r_eps)
 namespace coreneuron { } using namespace coreneuron; 
/*VERBATIM*/
/* ---------------------------------------------------------------- */
/* gammln - compute natural log of gamma function of xx */
static double
gammln(double xx)
{
    double x,tmp,ser;
    static double cof[6]={76.18009173,-86.50532033,24.01409822,
        -1.231739516,0.120858003e-2,-0.536382e-5};
    int j;
    x=xx-1.0;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.0;
    for (j=0;j<=5;j++) {
        x += 1.0;
        ser += cof[j]/x;
    }
    return -tmp+log(2.50662827465*ser);
}
 namespace coreneuron { 
double BnlDev ( _threadargsprotocomma_ double _lppr , double _lnnr ) {
   double _lBnlDev;
 
/*VERBATIM*/
        int j;
        double am,em,g,angle,p,bnl,sq,bt,y;
        double pc,plog,pclog,en,oldg;

        /* prepare to always ignore errors within this routine */

        p=(_lppr <= 0.5 ? _lppr : 1.0-_lppr);
        am=_lnnr*p;
        if (_lnnr < 25) {
            bnl=0.0;
            for (j=1;j<=_lnnr;j++)
                if (urand(_threadargs_) < p) bnl += 1.0;
        }
        else if (am < 1.0) {
            g=exp(-am);
            bt=1.0;
            for (j=0;j<=_lnnr;j++) {
                bt *= urand(_threadargs_);
                if (bt < g) break;
            }
            bnl=(j <= _lnnr ? j : _lnnr);
        }
        else {
            {
                en=_lnnr;
                oldg=gammln(en+1.0);
            }
            {
                pc=1.0-p;
                plog=log(p);
                pclog=log(pc);
            }
            sq=sqrt(2.0*am*pc);
            do {
                do {
                    angle=PI*urand(_threadargs_);
                    angle=PI*urand(_threadargs_);
                    y=tan(angle);
                    em=sq*y+am;
                } while (em < 0.0 || em >= (en+1.0));
                em=floor(em);
                    bt=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0) -
                    gammln(en-em+1.0)+em*plog+(en-em)*pclog);
            } while (urand(_threadargs_) > bt);
            bnl=em;
        }
        if (p != _lppr) bnl=_lnnr-bnl;

        /* recover error if changed during this routine, thus ignoring
            any errors during this routine */


        return bnl;

 _lBnlDev = bnl ;
   
return _lBnlDev;
 }
 
#if 0 /*BBCORE*/
 
static void _hoc_BnlDev(void) {
  double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  BnlDev ( _threadargs_, *getarg(1) , *getarg(2) ;
 hoc_retpushx(_r);
}
 
#endif /*BBCORE*/
 
double bbsavestate ( _threadargsproto_ ) {
   double _lbbsavestate;
 _lbbsavestate = 0.0 ;
   
/*VERBATIM*/
 #ifndef CORENEURON_BUILD
        // TODO: since N0,N1 are no longer state variables, they will need to be written using this callback
        //  provided that it is the version that supports multivalue writing
        /* first arg is direction (-1 get info, 0 save, 1 restore), second is value*/
        double *xdir, *xval, *hoc_pgetarg();
        long nrn_get_random_sequence(void* r);
        void nrn_set_random_sequence(void* r, int val);
        xdir = hoc_pgetarg(1);
        xval = hoc_pgetarg(2);
        if (_p_rng) {
                // tell how many items need saving
                if (*xdir == -1.) {
                    if( usingR123 ) {
                        *xdir = 2.0;
                    } else {
                        *xdir = 1.0;
                    }
                    return 0.0;
                }
                else if (*xdir == 0.) {
                    if( usingR123 ) {
                        uint32_t seq;
                        unsigned char which;
                        nrnran123_getseq( (nrnran123_State*)_p_rng, &seq, &which );
                        xval[0] = (double) seq;
                        xval[1] = (double) which;
                    } else {
                        xval[0] = (double)nrn_get_random_sequence(_p_rng);
                    }
                } else{
                    if( usingR123 ) {
                        nrnran123_setseq( (nrnran123_State*)_p_rng, (uint32_t)xval[0], (char)xval[1] );
                    } else {
                        nrn_set_random_sequence(_p_rng, (long)(xval[0]));
                    }
                }
        }

        // TODO: check for random123 and get the seq values
#endif
 
return _lbbsavestate;
 }
 
#if 0 /*BBCORE*/
 
static void _hoc_bbsavestate(void) {
  double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  bbsavestate ( _threadargs_ ;
 hoc_retpushx(_r);
}
 
#endif /*BBCORE*/
 static void _update_ion_pointer(Datum* _ppvar) {
 }

static inline void initmodel(_threadargsproto_) {
  int _i; double _save;{
  l = l0;
  n = n0;
 {
   
/*VERBATIM*/
    if (cvode_active_ && !deterministic) {
        hoc_execerror("StochKv2 with deterministic=0", "cannot be used with cvode");
    }

    if( usingR123 ) {
        nrnran123_setseq((nrnran123_State*)_p_rng, 0, 0);
    }
 eta = ( gkbar / gamma ) * ( 10000.0 ) ;
   trates ( _threadargscomma_ v ) ;
   n = ninf ;
   l = linf ;
   scale_dens = gamma / area ;
   N = floor ( eta * area + 0.5 ) ;
   N1L1 = n * l * N ;
   N1L0 = n * ( 1.0 - l ) * N ;
   N0L1 = ( 1.0 - n ) * l * N ;
   if (  ! deterministic ) {
     N1L1 = floor ( N1L1 + 0.5 ) ;
     N1L0 = floor ( N1L0 + 0.5 ) ;
     N0L1 = floor ( N0L1 + 0.5 ) ;
     }
   N0L0 = N - N1L1 - N1L0 - N0L1 ;
   n0l0_n1l0 = 0.0 ;
   n0l0_n0l1 = 0.0 ;
   n1l0_n1l1 = 0.0 ;
   n1l0_n0l0 = 0.0 ;
   n0l1_n1l1 = 0.0 ;
   n0l1_n0l0 = 0.0 ;
   n1l1_n0l1 = 0.0 ;
   n1l1_n1l0 = 0.0 ;
   }
 
}
}

void nrn_init(NrnThread* _nt, Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; ThreadDatum* _thread;
double _v, v; int* _ni; int _iml, _cntml_padded, _cntml_actual;
    _ni = _ml->_nodeindices;
_cntml_actual = _ml->_nodecount;
_cntml_padded = _ml->_nodecount_padded;
_thread = _ml->_thread;

#if defined(PG_ACC_BUGS)
#if defined(celsius)
#undef celsius;
_celsius_ = celsius;
#pragma acc update device (_celsius_) if(_nt->compute_gpu)
#define celsius _celsius_
#endif
#endif
_ACC_GLOBALS_UPDATE_
double * _nt_data = _nt->_data;
double * _vec_v = _nt->_actual_v;
int stream_id = _nt->stream_id;
  if (_nrn_skip_initmodel == 0) {
#if LAYOUT == 1 /*AoS*/
for (_iml = 0; _iml < _cntml_actual; ++_iml) {
 _p = _ml->_data + _iml*_psize; _ppvar = _ml->_pdata + _iml*_ppsize;
#elif LAYOUT == 0 /*SoA*/
 _p = _ml->_data; _ppvar = _ml->_pdata;
/* insert compiler dependent ivdep like pragma */
_PRAGMA_FOR_VECTOR_LOOP_
_PRAGMA_FOR_INIT_ACC_LOOP_
for (_iml = 0; _iml < _cntml_actual; ++_iml) {
#else /* LAYOUT > 1 */ /*AoSoA*/
#error AoSoA not implemented.
for (;;) { /* help clang-format properly indent */
#endif
    int _nd_idx = _ni[_iml];

#if 0
 _check_trates(_threadargs_);
#endif
    _v = _vec_v[_nd_idx];
    _PRCELLSTATE_V
 v = _v;
 _PRCELLSTATE_V
  ek = _ion_ek;
 initmodel(_threadargs_);
 }
  }
}

static double _nrn_current(_threadargsproto_, double _v){double _current=0.;v=_v;{ {
   gk = ( strap ( _threadargscomma_ N1L1 ) * scale_dens ) * ( 0.0001 ) ;
   ik = gk * ( v - ek ) ;
   }
 _current += ik;

} return _current;
}

#if defined(ENABLE_CUDA_INTERFACE) && defined(_OPENACC)
  void nrn_state_launcher(NrnThread*, Memb_list*, int, int);
  void nrn_jacob_launcher(NrnThread*, Memb_list*, int, int);
  void nrn_cur_launcher(NrnThread*, Memb_list*, int, int);
#endif


void nrn_cur(NrnThread* _nt, Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; ThreadDatum* _thread;
int* _ni; double _rhs, _g, _v, v; int _iml, _cntml_padded, _cntml_actual;
    _ni = _ml->_nodeindices;
_cntml_actual = _ml->_nodecount;
_cntml_padded = _ml->_nodecount_padded;
_thread = _ml->_thread;
double * _vec_rhs = _nt->_actual_rhs;
double * _vec_d = _nt->_actual_d;

#if defined(ENABLE_CUDA_INTERFACE) && defined(_OPENACC) && !defined(DISABLE_OPENACC)
  NrnThread* d_nt = acc_deviceptr(_nt);
  Memb_list* d_ml = acc_deviceptr(_ml);
  nrn_cur_launcher(d_nt, d_ml, _type, _cntml_actual);
  return;
#endif

double * _nt_data = _nt->_data;
double * _vec_v = _nt->_actual_v;
int stream_id = _nt->stream_id;
#if LAYOUT == 1 /*AoS*/
for (_iml = 0; _iml < _cntml_actual; ++_iml) {
 _p = _ml->_data + _iml*_psize; _ppvar = _ml->_pdata + _iml*_ppsize;
#elif LAYOUT == 0 /*SoA*/
 _p = _ml->_data; _ppvar = _ml->_pdata;
/* insert compiler dependent ivdep like pragma */
_PRAGMA_FOR_VECTOR_LOOP_
_PRAGMA_FOR_CUR_ACC_LOOP_
for (_iml = 0; _iml < _cntml_actual; ++_iml) {
#else /* LAYOUT > 1 */ /*AoSoA*/
#error AoSoA not implemented.
for (;;) { /* help clang-format properly indent */
#endif
    int _nd_idx = _ni[_iml];
    _v = _vec_v[_nd_idx];
    _PRCELLSTATE_V
  ek = _ion_ek;
 _g = _nrn_current(_threadargs_, _v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_threadargs_, _v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
 _PRCELLSTATE_G
	_vec_rhs[_nd_idx] -= _rhs;
	_vec_d[_nd_idx] += _g;
 
}
 
}

void nrn_state(NrnThread* _nt, Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; ThreadDatum* _thread;
double v, _v = 0.0; int* _ni; int _iml, _cntml_padded, _cntml_actual;
    _ni = _ml->_nodeindices;
_cntml_actual = _ml->_nodecount;
_cntml_padded = _ml->_nodecount_padded;
_thread = _ml->_thread;

#if defined(ENABLE_CUDA_INTERFACE) && defined(_OPENACC) && !defined(DISABLE_OPENACC)
  NrnThread* d_nt = acc_deviceptr(_nt);
  Memb_list* d_ml = acc_deviceptr(_ml);
  nrn_state_launcher(d_nt, d_ml, _type, _cntml_actual);
  return;
#endif

double * _nt_data = _nt->_data;
double * _vec_v = _nt->_actual_v;
int stream_id = _nt->stream_id;
#if LAYOUT == 1 /*AoS*/
for (_iml = 0; _iml < _cntml_actual; ++_iml) {
 _p = _ml->_data + _iml*_psize; _ppvar = _ml->_pdata + _iml*_ppsize;
#elif LAYOUT == 0 /*SoA*/
 _p = _ml->_data; _ppvar = _ml->_pdata;
/* insert compiler dependent ivdep like pragma */
_PRAGMA_FOR_VECTOR_LOOP_
_PRAGMA_FOR_STATE_ACC_LOOP_
for (_iml = 0; _iml < _cntml_actual; ++_iml) {
#else /* LAYOUT > 1 */ /*AoSoA*/
#error AoSoA not implemented.
for (;;) { /* help clang-format properly indent */
#endif
    int _nd_idx = _ni[_iml];
    _v = _vec_v[_nd_idx];
    _PRCELLSTATE_V
 v=_v;
{
  ek = _ion_ek;
 {   states(_threadargs_);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
 int _cntml_actual=1;
 int _cntml_padded=1;
 int _iml=0;
  if (!_first) return;
 
 _slist1 = (int*)malloc(sizeof(int)*2);
 _dlist1 = (int*)malloc(sizeof(int)*2);
 _slist1[0] = &(l) - _p;  _dlist1[0] = &(Dl) - _p;
 _slist1[1] = &(n) - _p;  _dlist1[1] = &(Dn) - _p;
 #pragma acc enter data copyin(_slist1[0:2])
 #pragma acc enter data copyin(_dlist1[0:2])

   _t_ntau = makevector(200*sizeof(double));
   _t_ltau = makevector(200*sizeof(double));
   _t_ninf = makevector(200*sizeof(double));
   _t_linf = makevector(200*sizeof(double));
   _t_al = makevector(200*sizeof(double));
   _t_bl = makevector(200*sizeof(double));
   _t_an = makevector(200*sizeof(double));
   _t_bn = makevector(200*sizeof(double));
_first = 0;
}
} // namespace coreneuron_lib
