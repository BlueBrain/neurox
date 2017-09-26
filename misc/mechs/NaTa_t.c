/* Created by Language version: 6.2.0 */
/* VECTORIZED */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "coreneuron/mech/cfile/scoplib.h"
#undef PI
 
#include "coreneuron/nrnoc/md1redef.h"
#include "coreneuron/coreneuron.h"

#if defined(_OPENACC) && !defined(DISABLE_OPENACC)
#include "coreneuron/nrniv/nrn_acc_manager.h"
#endif

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
 
#define nrn_init _nrn_init__NaTa_t
#define nrn_cur _nrn_cur__NaTa_t
#define nrn_cur_parallel _nrn_cur_parallel__NaTa_t
#define _nrn_current _nrn_current__NaTa_t
#define nrn_jacob _nrn_jacob__NaTa_t
#define nrn_state _nrn_state__NaTa_t
#define initmodel initmodel__NaTa_t
#define _net_receive _net_receive__NaTa_t
#define nrn_state_launcher nrn_state_NaTa_t_launcher
#define nrn_cur_launcher nrn_cur_NaTa_t_launcher
#define nrn_jacob_launcher nrn_jacob_NaTa_t_launcher 
#define rates rates_NaTa_t 
#define states states_NaTa_t 
#define _ode_matsol1 _nrn_ode_matsol1__NaTa_t
#define _ode_spec1 _nrn_ode_spec1__NaTa_t
 
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
#define gNaTa_tbar _p[0*_STRIDE]
#define ina _p[1*_STRIDE]
#define gNaTa_t _p[2*_STRIDE]
#define m _p[3*_STRIDE]
#define h _p[4*_STRIDE]
#define ttxo _p[5*_STRIDE]
#define ttxi _p[6*_STRIDE]
#define ena _p[7*_STRIDE]
#define mInf _p[8*_STRIDE]
#define mTau _p[9*_STRIDE]
#define mAlpha _p[10*_STRIDE]
#define mBeta _p[11*_STRIDE]
#define hInf _p[12*_STRIDE]
#define hTau _p[13*_STRIDE]
#define hAlpha _p[14*_STRIDE]
#define hBeta _p[15*_STRIDE]
#define Dm _p[16*_STRIDE]
#define Dh _p[17*_STRIDE]
#define _v_unused _p[18*_STRIDE]
#define _g_unused _p[19*_STRIDE]
 
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
#define _ion_ena		_nt_data[_ppvar[0*_STRIDE]]
#define _ion_ina	_nt_data[_ppvar[1*_STRIDE]]
#define _ion_dinadv	_nt_data[_ppvar[2*_STRIDE]]
#define _ion_ttxo		_nt_data[_ppvar[3*_STRIDE]]
#define _ion_ttxi		_nt_data[_ppvar[4*_STRIDE]]
 
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
 static int hoc_nrnpointerindex =  -1;
 static ThreadDatum* _extcall_thread;
 /* external NEURON variables */
 
#if 0 /*BBCORE*/
 /* declaration of user functions */
 static void _hoc_rates(void);
 
#endif /*BBCORE*/
 static int _mechtype;
 extern int nrn_get_mechtype();
extern void hoc_register_prop_size(int, int, int);
extern Memb_func* memb_func;
 
#if 0 /*BBCORE*/
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_NaTa_t", _hoc_setdata,
 "rates_NaTa_t", _hoc_rates,
 0, 0
};
 
#endif /*BBCORE*/
 /* declare global and static user variables */
 
static void _acc_globals_update() {
 }
 
#if 0 /*BBCORE*/
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gNaTa_tbar_NaTa_t", "S/cm2",
 "ina_NaTa_t", "mA/cm2",
 "gNaTa_t_NaTa_t", "S/cm2",
 0,0
};
 
#endif /*BBCORE*/
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
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
"NaTa_t",
 "gNaTa_tbar_NaTa_t",
 0,
 "ina_NaTa_t",
 "gNaTa_t_NaTa_t",
 0,
 "m_NaTa_t",
 "h_NaTa_t",
 0,
 0};
 static int _na_type;
 static int _ttx_type;
 
void _nrn_ode_state_vars__NaTa_t(short * count, short** var_offsets, short ** dv_offsets)
{
    *count = 2;
    (*var_offsets) = (short*) malloc(sizeof(short)* *count);
    (*dv_offsets) = (short*) malloc(sizeof(short)* *count);
    (*var_offsets)[0] = 3;
    (*var_offsets)[1] = 4;
    (*dv_offsets)[0] = 16;
    (*dv_offsets)[1] = 17;
}

static void nrn_alloc(double* _p, Datum* _ppvar, int _type) {
 
#if 0 /*BBCORE*/
 	/*initialize range parameters*/
 	gNaTa_tbar = 1e-05;
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 prop_ion = need_memb(_ttx_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[3]._pval = &prop_ion->param[2]; /* ttxo */
 	_ppvar[4]._pval = &prop_ion->param[1]; /* ttxi */
 
#endif /* BBCORE */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 
#define _psize 20
#define _ppsize 5
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*f)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(_threadargsproto_, int));
extern void _cvode_abstol( Symbol**, double*, int);

 void _NaTa_t_reg() {
	int _vectorized = 1;
  _initlists();
 _mechtype = nrn_get_mechtype(_mechanism[1]);
 if (_mechtype == -1) return;
 _nrn_layout_reg(_mechtype, LAYOUT);
 _na_type = nrn_get_mechtype("na_ion"); _ttx_type = nrn_get_mechtype("ttx_ion"); 
#if 0 /*BBCORE*/
 	ion_reg("na", -10000.);
 	ion_reg("ttx", 1.0);
 	_na_sym = hoc_lookup("na_ion");
 	_ttx_sym = hoc_lookup("ttx_ion");
 
#endif /*BBCORE*/
 	register_mech(_mechanism, nrn_alloc,nrn_cur, NULL, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
  hoc_register_prop_size(_mechtype, _psize, _ppsize);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ttx_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "ttx_ion");
 	hoc_register_var(hoc_scdoub, hoc_vdoub, NULL);
 }
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsproto_);
 
int _ode_spec1(_threadargsproto_);
/*int _ode_matsol1(_threadargsproto_);*/
 
#define _slist1 _slist1_NaTa_t
int* _slist1;
#pragma acc declare create(_slist1)

#define _dlist1 _dlist1_NaTa_t
int* _dlist1;
#pragma acc declare create(_dlist1)
 static inline int states(_threadargsproto_);
 
/*CVODE*/
int _ode_spec1 (_threadargsproto_) {int _reset = 0; {
   if ( ttxi  == 0.015625  && ttxo > 1e-12 ) {
     mInf = 0.0 ;
     mTau = 1e-12 ;
     hInf = 1.0 ;
     hTau = 1e-12 ;
     }
   else {
     rates ( _threadargs_ ) ;
     }
   Dm = ( mInf - m ) / mTau ;
   Dh = ( hInf - h ) / hTau ;
   }
 return _reset;
}
int _ode_matsol1 (_threadargsproto_) {
 if ( ttxi  == 0.015625  && ttxo > 1e-12 ) {
   mInf = 0.0 ;
   mTau = 1e-12 ;
   hInf = 1.0 ;
   hTau = 1e-12 ;
   }
 else {
   rates ( _threadargs_ ) ;
   }
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mTau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / hTau )) ;
 return 0;
}
 /*END CVODE*/
 static int states (_threadargsproto_) { {
   if ( ttxi  == 0.015625  && ttxo > 1e-12 ) {
     mInf = 0.0 ;
     mTau = 1e-12 ;
     hInf = 1.0 ;
     hTau = 1e-12 ;
     }
   else {
     rates ( _threadargs_ ) ;
     }
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / mTau)))*(- ( ( ( mInf ) ) / mTau ) / ( ( ( ( - 1.0) ) ) / mTau ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / hTau)))*(- ( ( ( hInf ) ) / hTau ) / ( ( ( ( - 1.0) ) ) / hTau ) - h) ;
   }
  return 0;
}
 
static int  rates ( _threadargsproto_ ) {
   double _lqt ;
 _lqt = pow( 2.3 , ( ( 34.0 - 21.0 ) / 10.0 ) ) ;
    if ( v  == - 38.0 ) {
     v = v + 0.0001 ;
     }
   mAlpha = ( 0.182 * ( v - - 38.0 ) ) / ( 1.0 - ( exp ( - ( v - - 38.0 ) / 6.0 ) ) ) ;
   mBeta = ( 0.124 * ( - v - 38.0 ) ) / ( 1.0 - ( exp ( - ( - v - 38.0 ) / 6.0 ) ) ) ;
   mTau = ( 1.0 / ( mAlpha + mBeta ) ) / _lqt ;
   mInf = mAlpha / ( mAlpha + mBeta ) ;
   if ( v  == - 66.0 ) {
     v = v + 0.0001 ;
     }
   hAlpha = ( - 0.015 * ( v - - 66.0 ) ) / ( 1.0 - ( exp ( ( v - - 66.0 ) / 6.0 ) ) ) ;
   hBeta = ( - 0.015 * ( - v - 66.0 ) ) / ( 1.0 - ( exp ( ( - v - 66.0 ) / 6.0 ) ) ) ;
   hTau = ( 1.0 / ( hAlpha + hBeta ) ) / _lqt ;
   hInf = hAlpha / ( hAlpha + hBeta ) ;
     return 0; }
 
#if 0 /*BBCORE*/
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rates ( _threadargs_ ;
 hoc_retpushx(_r);
}
 
#endif /*BBCORE*/
 static void _update_ion_pointer(Datum* _ppvar) {
 }

static void initmodel(_threadargsproto_) {
  int _i; double _save;{
  h = h0;
  m = m0;
 {
   if ( ttxi  == 0.015625  && ttxo > 1e-12 ) {
     mInf = 0.0 ;
     mTau = 1e-12 ;
     hInf = 1.0 ;
     hTau = 1e-12 ;
     }
   else {
     rates ( _threadargs_ ) ;
     }
   m = mInf ;
   h = hInf ;
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
  ena = _ion_ena;
  ttxo = _ion_ttxo;
  ttxi = _ion_ttxi;
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
   gNaTa_t = gNaTa_tbar * m * m * m * h ;
   ina = gNaTa_t * ( v - ena ) ;
   }
 _current += ina;

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
                        const mod_acc_f_t acc_rhs_d, const mod_acc_f_t acc_i_didv, void *args)
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
  ena = _ion_ena;
  ttxo = _ion_ttxo;
  ttxi = _ion_ttxi;
 _g = _nrn_current(_threadargs_, _v + .001);
  double _dina;
  _dina = ina;
 _rhs = _nrn_current(_threadargs_, _v);
 _g = (_g - _rhs)/.001;
 if (acc_i_didv)
 {
   _vec_shadow_didv[_iml] = + (_dina - ina)/.001;
   _vec_shadow_i[_iml] = +ina;
 }
 else
 {
   _ion_dinadv += (_dina - ina)/.001 ;
   _ion_ina += ina ;
 }
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
  ena = _ion_ena;
  ttxo = _ion_ttxo;
  ttxi = _ion_ttxi;
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
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
 #pragma acc enter data copyin(_slist1[0:2])
 #pragma acc enter data copyin(_dlist1[0:2])

_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif
