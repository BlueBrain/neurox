/* Created by Language version: 6.2.0 */
/* VECTORIZED */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "coreneuron/mech/cfile/scoplib.h"
#undef PI
 
#include "coreneuron/nrnoc/md1redef.h"
#include "coreneuron/nrnconf.h"
#include "coreneuron/nrnoc/multicore.h"

#if defined(_OPENACC) && !defined(DISABLE_OPENACC)
#include "coreneuron/nrniv/nrn_acc_manager.h"

#endif
#include "coreneuron/utils/randoms/nrnran123.h"

#include "coreneuron/nrnoc/md2redef.h"
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#if !defined(DISABLE_HOC_EXP)
#undef exp
#define exp hoc_Exp
#endif
extern double hoc_Exp(double);
#endif
 
#define _thread_present_ /**/ , _thread[0:1] , _slist1[0:1], _dlist1[0:1] 
 
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
 
#define nrn_init _nrn_init__StochKv
#define nrn_cur _nrn_cur__StochKv
#define nrn_cur_parallel _nrn_cur_parallel__StochKv
#define _nrn_current _nrn_current__StochKv
#define nrn_jacob _nrn_jacob__StochKv
#define nrn_state _nrn_state__StochKv
#define initmodel initmodel__StochKv
#define _net_receive _net_receive__StochKv
#define nrn_state_launcher nrn_state_StochKv_launcher
#define nrn_cur_launcher nrn_cur_StochKv_launcher
#define nrn_jacob_launcher nrn_jacob_StochKv_launcher 
#define ChkProb ChkProb_StochKv 
#define setRNG setRNG_StochKv 
#define states states_StochKv 
#define trates trates_StochKv 
 
#define _threadargscomma_ _iml, _cntml_padded, _p, _ppvar, _thread, _nt, v,
#define _threadargsprotocomma_ int _iml, int _cntml_padded, double* _p, Datum* _ppvar, ThreadDatum* _thread, _NrnThread* _nt, double v,
#define _threadargs_ _iml, _cntml_padded, _p, _ppvar, _thread, _nt, v
#define _threadargsproto_ int _iml, int _cntml_padded, double* _p, Datum* _ppvar, ThreadDatum* _thread, _NrnThread* _nt, double v
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gamma _p[0*_STRIDE]
#define eta _p[1*_STRIDE]
#define gkbar _p[2*_STRIDE]
#define ik _p[3*_STRIDE]
#define gk _p[4*_STRIDE]
#define N _p[5*_STRIDE]
#define N0 _p[6*_STRIDE]
#define N1 _p[7*_STRIDE]
#define n0_n1 _p[8*_STRIDE]
#define n1_n0 _p[9*_STRIDE]
#define n _p[10*_STRIDE]
#define ek _p[11*_STRIDE]
#define scale_dens _p[12*_STRIDE]
#define n0_n1_new _p[13*_STRIDE]
#define usingR123 _p[14*_STRIDE]
#define Dn _p[15*_STRIDE]
#define _v_unused _p[16*_STRIDE]
#define _g_unused _p[17*_STRIDE]
 
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
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  3;
 static ThreadDatum* _extcall_thread;
 /* external NEURON variables */
 extern double celsius;
 #if defined(PG_ACC_BUGS)
#define _celsius_ _celsius__StochKv
double _celsius_;
#pragma acc declare copyin(_celsius_)
#define celsius _celsius_
#endif
 
#if 0 /*BBCORE*/
 /* declaration of user functions */
 static void _hoc_BnlDev(void);
 static void _hoc_ChkProb(void);
 static void _hoc_SigmoidRate(void);
 static void _hoc_bbsavestate(void);
 static void _hoc_brand(void);
 static void _hoc_setRNG(void);
 static void _hoc_strap(void);
 static void _hoc_trates(void);
 static void _hoc_urand(void);
 
#endif /*BBCORE*/
 static int _mechtype;
 extern int nrn_get_mechtype();
extern void hoc_register_prop_size(int, int, int);
extern Memb_func* memb_func;
 
#if 0 /*BBCORE*/
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_StochKv", _hoc_setdata,
 "BnlDev_StochKv", _hoc_BnlDev,
 "ChkProb_StochKv", _hoc_ChkProb,
 "SigmoidRate_StochKv", _hoc_SigmoidRate,
 "bbsavestate_StochKv", _hoc_bbsavestate,
 "brand_StochKv", _hoc_brand,
 "setRNG_StochKv", _hoc_setRNG,
 "strap_StochKv", _hoc_strap,
 "trates_StochKv", _hoc_trates,
 "urand_StochKv", _hoc_urand,
 0, 0
};
 
#endif /*BBCORE*/
#define BnlDev BnlDev_StochKv
#define SigmoidRate SigmoidRate_StochKv
#define bbsavestate bbsavestate_StochKv
#define brand brand_StochKv
#define strap strap_StochKv
#define urand urand_StochKv
 inline double BnlDev( _threadargsprotocomma_ double , double );
 inline double SigmoidRate( _threadargsprotocomma_ double , double , double , double );
 inline double bbsavestate( _threadargsproto_ );
 inline double brand( _threadargsprotocomma_ double , double );
 inline double strap( _threadargsprotocomma_ double );
 inline double urand( _threadargsproto_ );
 /* declare global and static user variables */
 static int _thread1data_inuse = 0;
static double _thread1data[7];
#define _gth 0
#define P_b_StochKv _thread1data[0]
#define P_b _thread[_gth]._pval[0]
#define P_a_StochKv _thread1data[1]
#define P_a _thread[_gth]._pval[1]
#define Rb Rb_StochKv
 double Rb = 0.002;
 #pragma acc declare copyin (Rb)
#define Ra Ra_StochKv
 double Ra = 0.02;
 #pragma acc declare copyin (Ra)
#define a_StochKv _thread1data[2]
#define a _thread[_gth]._pval[2]
#define b_StochKv _thread1data[3]
#define b _thread[_gth]._pval[3]
#define deterministic deterministic_StochKv
 double deterministic = 0;
 #pragma acc declare copyin (deterministic)
#define ntau_StochKv _thread1data[4]
#define ntau _thread[_gth]._pval[4]
#define ninf_StochKv _thread1data[5]
#define ninf _thread[_gth]._pval[5]
#define qa qa_StochKv
 double qa = 9;
 #pragma acc declare copyin (qa)
#define q10 q10_StochKv
 double q10 = 2.3;
 #pragma acc declare copyin (q10)
#define tha tha_StochKv
 double tha = -40;
 #pragma acc declare copyin (tha)
#define tadj_StochKv _thread1data[6]
#define tadj _thread[_gth]._pval[6]
#define temp temp_StochKv
 double temp = 23;
 #pragma acc declare copyin (temp)
#define vmax vmax_StochKv
 double vmax = 100;
 #pragma acc declare copyin (vmax)
#define vmin vmin_StochKv
 double vmin = -120;
 #pragma acc declare copyin (vmin)
 
static void _acc_globals_update() {
 #pragma acc update device (Rb) if(nrn_threads->compute_gpu)
 #pragma acc update device (Ra) if(nrn_threads->compute_gpu)
 #pragma acc update device (deterministic) if(nrn_threads->compute_gpu)
 #pragma acc update device (qa) if(nrn_threads->compute_gpu)
 #pragma acc update device (q10) if(nrn_threads->compute_gpu)
 #pragma acc update device (tha) if(nrn_threads->compute_gpu)
 #pragma acc update device (temp) if(nrn_threads->compute_gpu)
 #pragma acc update device (vmax) if(nrn_threads->compute_gpu)
 #pragma acc update device (vmin) if(nrn_threads->compute_gpu)
 }
 
#if 0 /*BBCORE*/
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tha_StochKv", "mV",
 "Ra_StochKv", "/ms",
 "Rb_StochKv", "/ms",
 "temp_StochKv", "degC",
 "vmin_StochKv", "mV",
 "vmax_StochKv", "mV",
 "a_StochKv", "/ms",
 "b_StochKv", "/ms",
 "ntau_StochKv", "ms",
 "gamma_StochKv", "pS",
 "eta_StochKv", "1/um2",
 "gkbar_StochKv", "S/cm2",
 "ik_StochKv", "mA/cm2",
 "gk_StochKv", "S/cm2",
 0,0
};
 
#endif /*BBCORE*/
 static double delta_t = 1;
 static double n0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "tha_StochKv", &tha_StochKv,
 "qa_StochKv", &qa_StochKv,
 "Ra_StochKv", &Ra_StochKv,
 "Rb_StochKv", &Rb_StochKv,
 "temp_StochKv", &temp_StochKv,
 "q10_StochKv", &q10_StochKv,
 "deterministic_StochKv", &deterministic_StochKv,
 "vmin_StochKv", &vmin_StochKv,
 "vmax_StochKv", &vmax_StochKv,
 "a_StochKv", &a_StochKv,
 "b_StochKv", &b_StochKv,
 "ninf_StochKv", &ninf_StochKv,
 "ntau_StochKv", &ntau_StochKv,
 "tadj_StochKv", &tadj_StochKv,
 "P_a_StochKv", &P_a_StochKv,
 "P_b_StochKv", &P_b_StochKv,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(double*, Datum*, int);
void nrn_init(_NrnThread*, _Memb_list*, int);
void nrn_state(_NrnThread*, _Memb_list*, int);
 void nrn_cur(_NrnThread*, _Memb_list*, int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "6.2.0",
"StochKv",
 "gamma_StochKv",
 "eta_StochKv",
 "gkbar_StochKv",
 0,
 "ik_StochKv",
 "gk_StochKv",
 "N_StochKv",
 "N0_StochKv",
 "N1_StochKv",
 "n0_n1_StochKv",
 "n1_n0_StochKv",
 0,
 "n_StochKv",
 0,
 "rng_StochKv",
 0};
 static int _k_type;
 
static void nrn_alloc(double* _p, Datum* _ppvar, int _type) {
 
#if 0 /*BBCORE*/
 	/*initialize range parameters*/
 	gamma = 30;
 	eta = 0;
 	gkbar = 0.75;
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
#endif /* BBCORE */
 
}
 static void _initlists();
 static void _thread_mem_init(ThreadDatum*);
 static void _thread_cleanup(ThreadDatum*);
 static void _update_ion_pointer(Datum*);
 
#define _psize 18
#define _ppsize 5
 static void bbcore_read(double *, int*, int*, int*, _threadargsproto_);
 extern void hoc_reg_bbcore_read(int, void(*)(double *, int*, int*, int*, _threadargsproto_));
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*f)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(_threadargsproto_, int));
extern void _cvode_abstol( Symbol**, double*, int);

 void _StochKv_reg() {
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
 	register_mech(_mechanism, nrn_alloc,nrn_cur, NULL, nrn_state, nrn_init, hoc_nrnpointerindex, 2);
  _extcall_thread = (ThreadDatum*)ecalloc(1, sizeof(ThreadDatum));
  _thread_mem_init(_extcall_thread);
  _thread1data_inuse = 0;
     _nrn_thread_reg1(_mechtype, _thread_mem_init);
     _nrn_thread_reg0(_mechtype, _thread_cleanup);
   hoc_reg_bbcore_read(_mechtype, bbcore_read);
  hoc_register_prop_size(_mechtype, _psize, _ppsize);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "bbcorepointer");
  hoc_register_dparam_semantics(_mechtype, 4, "area");
 	hoc_register_var(hoc_scdoub, hoc_vdoub, NULL);
 }
static char *modelname = "skm95.mod  ";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int ChkProb(_threadargsprotocomma_ double);
static int setRNG(_threadargsproto_);
static int trates(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 
#define _slist1 _slist1_StochKv
int* _slist1;
#pragma acc declare create(_slist1)

#define _dlist1 _dlist1_StochKv
int* _dlist1;
#pragma acc declare create(_dlist1)
 static inline int states(_threadargsproto_);
 
/*VERBATIM*/
#include "nrnran123.h"
extern int cvode_active_;

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#if !defined(CORENEURON_BUILD)
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
#endif

 
/*CVODE*/
 static int _ode_spec1 (_threadargsproto_) {int _reset = 0; {
   trates ( _threadargscomma_ v ) ;
   Dn = a - ( a + b ) * n ;
   if ( deterministic  || dt > 1.0 ) {
     N1 = n * N ;
     }
   else {
     P_a = strap ( _threadargscomma_ a * dt ) ;
     P_b = strap ( _threadargscomma_ b * dt ) ;
     ChkProb ( _threadargscomma_ P_a ) ;
     ChkProb ( _threadargscomma_ P_b ) ;
     n0_n1 = BnlDev ( _threadargscomma_ P_a , N0 ) ;
     n1_n0 = BnlDev ( _threadargscomma_ P_b , N1 ) ;
     N0 = strap ( _threadargscomma_ N0 - n0_n1 + n1_n0 ) ;
     N1 = N - N0 ;
     }
   N0 = N - N1 ;
   }
 return _reset;
}
 static int _ode_matsol1 (_threadargsproto_) {
 trates ( _threadargscomma_ v ) ;
 Dn = Dn  / (1. - dt*( ( - (( a + b ))*(1.0) ) )) ;
 return 0;
}
 /*END CVODE*/
 static int states (_threadargsproto_) { {
   trates ( _threadargscomma_ v ) ;
    n = n + (1. - exp(dt*(( - (( a + b ))*(1.0) ))))*(- ( a ) / ( ( - (( a + b ))*(1.0)) ) - n) ;
   if ( deterministic  || dt > 1.0 ) {
     N1 = n * N ;
     }
   else {
     P_a = strap ( _threadargscomma_ a * dt ) ;
     P_b = strap ( _threadargscomma_ b * dt ) ;
     ChkProb ( _threadargscomma_ P_a ) ;
     ChkProb ( _threadargscomma_ P_b ) ;
     n0_n1 = BnlDev ( _threadargscomma_ P_a , N0 ) ;
     n1_n0 = BnlDev ( _threadargscomma_ P_b , N1 ) ;
     N0 = strap ( _threadargscomma_ N0 - n0_n1 + n1_n0 ) ;
     N1 = N - N0 ;
     }
   N0 = N - N1 ;
   }
  return 0;
}
 
static int  trates ( _threadargsprotocomma_ double _lv ) {
   tadj = pow( q10 , ( ( celsius - temp ) / ( 10.0 ) ) ) ;
   a = SigmoidRate ( _threadargscomma_ _lv , tha , Ra , qa ) ;
   a = a * tadj ;
   b = SigmoidRate ( _threadargscomma_ - _lv , - tha , Rb , qa ) ;
   b = b * tadj ;
   ntau = 1.0 / ( a + b ) ;
   ninf = a * ntau ;
    return 0; }
 
#if 0 /*BBCORE*/
 
static void _hoc_trates(void) {
  double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 trates ( _threadargs_, *getarg(1) ;
 hoc_retpushx(_r);
}
 
#endif /*BBCORE*/
 
double SigmoidRate ( _threadargsprotocomma_ double _lv , double _lth , double _la , double _lq ) {
   double _lSigmoidRate;
  if ( fabs ( _lv - _lth ) > 1e-6 ) {
     _lSigmoidRate = _la * ( _lv - _lth ) / ( 1.0 - exp ( - ( _lv - _lth ) / _lq ) ) ;
      }
   else {
     _lSigmoidRate = _la * _lq ;
     }
   
return _lSigmoidRate;
 }
 
#if 0 /*BBCORE*/
 
static void _hoc_SigmoidRate(void) {
  double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  SigmoidRate ( _threadargs_, *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) ;
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
   double* _p; Datum* _ppvar; ThreadDatum* _thread; _NrnThread* _nt;
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
    fprintf(stderr, "StochKv.mod:ChkProb: argument not a probability.\n");
 }
    return 0; }
 
#if 0 /*BBCORE*/
 
static void _hoc_ChkProb(void) {
  double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; _NrnThread* _nt;
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
#if !defined(NRNBBCORE) || !NRNBBCORE
    usingR123 = 0;
    if( ifarg(1) && hoc_is_double_arg(1) ) {
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
        uint32_t a2 = 0;
        
        if (*pv) {
            nrnran123_deletestream(*pv);
            *pv = (nrnran123_State*)0;
        } 
        if (ifarg(2)) {
            a2 = (uint32_t)*getarg(2);
        }
        *pv = nrnran123_newstream((uint32_t)*getarg(1), a2);
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
   double* _p; Datum* _ppvar; ThreadDatum* _thread; _NrnThread* _nt;
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
#if !defined(CORENEURON_BUILD)
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
   double* _p; Datum* _ppvar; ThreadDatum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  urand ( _threadargs_ ;
 hoc_retpushx(_r);
}
 
#endif /*BBCORE*/
 
/*VERBATIM*/
static void bbcore_write(double* x, int* d, int* xx, int* offset, _threadargsproto_) {
    if (d) {
        uint32_t* di = ((uint32_t*)d) + *offset;
      // temporary just enough to see how much space is being used
      if (!_p_rng) {
        di[0] = 0; di[1] = 0;
      }else{
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
        nrnran123_getids(*pv, di, di+1);
      }
//printf("StochKv.mod %p: bbcore_write offset=%d %d %d\n", _p, *offset, d?di[0]:-1, d?di[1]:-1);
    }
    *offset += 2;
}
static void bbcore_read(double* x, int* d, int* xx, int* offset, _threadargsproto_) {
    assert(!_p_rng);
    uint32_t* di = ((uint32_t*)d) + *offset;
        if (di[0] != 0 || di[1] != 0)
        {
      nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
      *pv = nrnran123_newstream(di[0], di[1]);
        }
//printf("StochKv.mod %p: bbcore_read offset=%d %d %d\n", _p, *offset, di[0], di[1]);
    *offset += 2;
}
 
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
   double* _p; Datum* _ppvar; ThreadDatum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  brand ( _threadargs_, *getarg(1) , *getarg(2) ;
 hoc_retpushx(_r);
}
 
#endif /*BBCORE*/
 
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
   double* _p; Datum* _ppvar; ThreadDatum* _thread; _NrnThread* _nt;
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
#ifdef ENABLE_SAVE_STATE
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
                        char which;
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
   double* _p; Datum* _ppvar; ThreadDatum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  bbsavestate ( _threadargs_ ;
 hoc_retpushx(_r);
}
 
#endif /*BBCORE*/
 
static void _thread_mem_init(ThreadDatum* _thread) {
  if (_thread1data_inuse) {_thread[_gth]._pval = (double*)ecalloc(7, sizeof(double));
 }else{
 _thread[_gth]._pval = _thread1data; _thread1data_inuse = 1;
 }
 }
 
static void _thread_cleanup(ThreadDatum* _thread) {
  if (_thread[_gth]._pval == _thread1data) {
   _thread1data_inuse = 0;
  }else{
   free((void*)_thread[_gth]._pval);
  }
 }
 static void _update_ion_pointer(Datum* _ppvar) {
 }

static void initmodel(_threadargsproto_) {
  int _i; double _save;{
  n = n0;
 {
   
/*VERBATIM*/
    if (_p_rng) {
        nrnran123_setseq((nrnran123_State*)_p_rng, 0, 0);
    }
    if (cvode_active_ && !deterministic) {
        hoc_execerror("StochKv with deterministic=0", "cannot be used with cvode");
    }
 eta = gkbar / gamma ;
   trates ( _threadargscomma_ v ) ;
   n = ninf ;
   scale_dens = gamma / area ;
   N = floor ( eta * area + 0.5 ) ;
   N1 = n * N ;
   if (  ! deterministic ) {
     N1 = floor ( N1 + 0.5 ) ;
     }
   N0 = N - N1 ;
   n0_n1 = 0.0 ;
   n1_n0 = 0.0 ;
   }
 
}
}

void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
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
    _v = _vec_v[_nd_idx];
    _PRCELLSTATE_V
 v = _v;
 _PRCELLSTATE_V
  ek = _ion_ek;
 initmodel(_threadargs_);

 //populate offsets arrays //(if parallel processing)
 if (_ml->_shadow_didv_offsets)
 {
   _ml->_shadow_i_offsets[_iml] = _ppvar[1*_STRIDE];
   _ml->_shadow_didv_offsets[_iml] = _ppvar[2*_STRIDE];
 }
 }
}

static double _nrn_current(_threadargsproto_, double _v){double _current=0.;v=_v;{ {
   gk = ( strap ( _threadargscomma_ N1 ) * scale_dens * tadj ) ;
   ik = 1e-4 * gk * ( v - ek ) ;
   }
 _current += ik;

} return _current;
}

#if defined(ENABLE_CUDA_INTERFACE) && defined(_OPENACC)
  void nrn_state_launcher(_NrnThread*, _Memb_list*, int, int);
  void nrn_jacob_launcher(_NrnThread*, _Memb_list*, int, int);
  void nrn_cur_launcher(_NrnThread*, _Memb_list*, int, int);
#endif

void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
    nrn_cur_parallel(_nt, _ml, _type, NULL, NULL, NULL);
}

void nrn_cur_parallel(_NrnThread* _nt, _Memb_list* _ml, int _type,
                        mod_acc_f_t acc_rhs_d, mod_acc_f_t acc_i_didv, void *args)
{
double* _p; Datum* _ppvar; ThreadDatum* _thread;
int* _ni; double _rhs, _g, _v, v; int _iml, _cntml_padded, _cntml_actual;
    _ni = _ml->_nodeindices;
_cntml_actual = _ml->_nodecount;
_cntml_padded = _ml->_nodecount_padded;
_thread = _ml->_thread;
double * _vec_rhs = _nt->_actual_rhs;
double * _vec_d =   _nt->_actual_d;
double * _vec_shadow_rhs = _ml->_shadow_rhs;
double * _vec_shadow_d = _ml->_shadow_d;
double * _vec_shadow_i = _ml->_shadow_i;
double * _vec_shadow_didv = _ml->_shadow_didv;


#if defined(ENABLE_CUDA_INTERFACE) && defined(_OPENACC) && !defined(DISABLE_OPENACC)
  _NrnThread* d_nt = acc_deviceptr(_nt);
  _Memb_list* d_ml = acc_deviceptr(_ml);
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
 if (acc_rhs_d)
 {
    _vec_shadow_rhs[_iml] = -_rhs;
    _vec_shadow_d[_iml] = +_g;
 }
 else
 {
    _vec_rhs[_nd_idx] -= _rhs;
    _vec_d[_nd_idx] += _g;
 }
}
//accumulation of individual contributions (for parallel executions)
if (acc_rhs_d)  (*acc_rhs_d) (_nt, _ml, _type, args);
if (acc_i_didv) (*acc_i_didv)(_nt, _ml, _type, args);
}

void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; ThreadDatum* _thread;
double v, _v = 0.0; int* _ni; int _iml, _cntml_padded, _cntml_actual;
    _ni = _ml->_nodeindices;
_cntml_actual = _ml->_nodecount;
_cntml_padded = _ml->_nodecount_padded;
_thread = _ml->_thread;

#if defined(ENABLE_CUDA_INTERFACE) && defined(_OPENACC) && !defined(DISABLE_OPENACC)
  _NrnThread* d_nt = acc_deviceptr(_nt);
  _Memb_list* d_ml = acc_deviceptr(_ml);
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
 
 _slist1 = (int*)malloc(sizeof(int)*1);
 _dlist1 = (int*)malloc(sizeof(int)*1);
 _slist1[0] = &(n) - _p;  _dlist1[0] = &(Dn) - _p;
 #pragma acc enter data copyin(_slist1[0:1])
 #pragma acc enter data copyin(_dlist1[0:1])

_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif
