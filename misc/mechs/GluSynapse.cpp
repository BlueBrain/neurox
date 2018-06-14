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
 
#define _thread_present_ /**/ , _slist1[0:10], _dlist1[0:10] 
 
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
 

#if !defined(NET_RECEIVE_BUFFERING)
#define NET_RECEIVE_BUFFERING 1
#endif
 
#define nrn_init _nrn_init__GluSynapse
#define nrn_cur _nrn_cur__GluSynapse
#define _nrn_current _nrn_current__GluSynapse
#define nrn_jacob _nrn_jacob__GluSynapse
#define nrn_state _nrn_state__GluSynapse
#define initmodel initmodel__GluSynapse
#define _net_receive _net_receive__GluSynapse
#define nrn_state_launcher nrn_state_GluSynapse_launcher
#define nrn_cur_launcher nrn_cur_GluSynapse_launcher
#define nrn_jacob_launcher nrn_jacob_GluSynapse_launcher 
#if NET_RECEIVE_BUFFERING
#define _net_buf_receive _net_buf_receive_GluSynapse
static void _net_buf_receive(NrnThread*);
#endif
 
#define _nrn_watch_check _nrn_watch_check__GluSynapse 
#define setRNG setRNG_GluSynapse 
#define state state_GluSynapse 
#define _ode_matsol1 _nrn_ode_matsol1__GluSynapse
#define _ode_spec1 _nrn_ode_spec1__GluSynapse

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
#define tau_r_AMPA _p[0*_STRIDE]
#define tau_d_AMPA _p[1*_STRIDE]
#define E_AMPA _p[2*_STRIDE]
#define gmax_AMPA _p[3*_STRIDE]
#define tau_r_NMDA _p[4*_STRIDE]
#define tau_d_NMDA _p[5*_STRIDE]
#define E_NMDA _p[6*_STRIDE]
#define Use _p[7*_STRIDE]
#define Dep _p[8*_STRIDE]
#define Fac _p[9*_STRIDE]
#define Nrrp _p[10*_STRIDE]
#define volume_CR _p[11*_STRIDE]
#define gca_bar_VDCC _p[12*_STRIDE]
#define ljp_VDCC _p[13*_STRIDE]
#define vhm_VDCC _p[14*_STRIDE]
#define km_VDCC _p[15*_STRIDE]
#define vhh_VDCC _p[16*_STRIDE]
#define kh_VDCC _p[17*_STRIDE]
#define tm_VDCC _p[18*_STRIDE]
#define th_VDCC _p[19*_STRIDE]
#define gamma_ca_CR _p[20*_STRIDE]
#define tau_ca_CR _p[21*_STRIDE]
#define min_ca_CR _p[22*_STRIDE]
#define cao_CR _p[23*_STRIDE]
#define tau_GB _p[24*_STRIDE]
#define theta_d_GB _p[25*_STRIDE]
#define theta_p_GB _p[26*_STRIDE]
#define gamma_d_GB _p[27*_STRIDE]
#define gamma_p_GB _p[28*_STRIDE]
#define rho_star_GB _p[29*_STRIDE]
#define rho0_GB _p[30*_STRIDE]
#define enable_GB _p[31*_STRIDE]
#define tau_Use_GB _p[32*_STRIDE]
#define Use_d_GB _p[33*_STRIDE]
#define Use_p_GB _p[34*_STRIDE]
#define NMDA_ratio _p[35*_STRIDE]
#define mg _p[36*_STRIDE]
#define synapseID _p[37*_STRIDE]
#define verbose _p[38*_STRIDE]
#define selected_for_report _p[39*_STRIDE]
#define factor_AMPA _p[40*_STRIDE]
#define g_AMPA _p[41*_STRIDE]
#define i_AMPA _p[42*_STRIDE]
#define factor_NMDA _p[43*_STRIDE]
#define gmax_NMDA _p[44*_STRIDE]
#define g_NMDA _p[45*_STRIDE]
#define i_NMDA _p[46*_STRIDE]
#define u _p[47*_STRIDE]
#define tsyn _p[48*_STRIDE]
#define Psurv _p[49*_STRIDE]
#define unoccupied _p[50*_STRIDE]
#define occupied _p[51*_STRIDE]
#define Pf_NMDA _p[52*_STRIDE]
#define ica_NMDA _p[53*_STRIDE]
#define area_CR _p[54*_STRIDE]
#define gca_VDCC _p[55*_STRIDE]
#define gca_bar_abs_VDCC _p[56*_STRIDE]
#define minf_VDCC _p[57*_STRIDE]
#define hinf_VDCC _p[58*_STRIDE]
#define mtau_VDCC _p[59*_STRIDE]
#define htau_VDCC _p[60*_STRIDE]
#define ica_VDCC _p[61*_STRIDE]
#define depress_GB _p[62*_STRIDE]
#define potentiate_GB _p[63*_STRIDE]
#define w _p[64*_STRIDE]
#define g _p[65*_STRIDE]
#define i _p[66*_STRIDE]
#define A_AMPA _p[67*_STRIDE]
#define B_AMPA _p[68*_STRIDE]
#define A_NMDA _p[69*_STRIDE]
#define B_NMDA _p[70*_STRIDE]
#define m_VDCC _p[71*_STRIDE]
#define h_VDCC _p[72*_STRIDE]
#define cai_CR _p[73*_STRIDE]
#define Rho_GB _p[74*_STRIDE]
#define Use_GB _p[75*_STRIDE]
#define effcai_GB _p[76*_STRIDE]
#define usingR123 _p[77*_STRIDE]
#define DA_AMPA _p[78*_STRIDE]
#define DB_AMPA _p[79*_STRIDE]
#define DA_NMDA _p[80*_STRIDE]
#define DB_NMDA _p[81*_STRIDE]
#define Dm_VDCC _p[82*_STRIDE]
#define Dh_VDCC _p[83*_STRIDE]
#define Dcai_CR _p[84*_STRIDE]
#define DRho_GB _p[85*_STRIDE]
#define DUse_GB _p[86*_STRIDE]
#define Deffcai_GB _p[87*_STRIDE]
#define _v_unused _p[88*_STRIDE]
#define _g_unused _p[89*_STRIDE]
#define _tsav _p[90*_STRIDE]
 
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
#define _nd_area  _nt_data[_ppvar[0*_STRIDE]]
#define _p_rng_rel	_nt->_vdata[_ppvar[2*_STRIDE]]
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 static int hoc_nrnpointerindex =  2;
 static ThreadDatum* _extcall_thread;
 #define FARADAY FARADAY_GluSynapse
 #define PI PI_GluSynapse
 #define R R_GluSynapse
 /* external NEURON variables */
 
#if 0 /*BBCORE*/
 /* declaration of user functions */
 static double _hoc_bbsavestate();
 static double _hoc_nernst();
 static double _hoc_setRNG();
 static double _hoc_toggleVerbose();
 static double _hoc_toggleLTPlasticity();
 static double _hoc_urand();
 
#endif /*BBCORE*/
 
#define _mechtype _mechtype_GluSynapse
int _mechtype;
#pragma acc declare copyin (_mechtype)
 static int _pointtype;
 
#if 0 /*BBCORE*/
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt();
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 
#endif /*BBCORE*/
 
#if 0 /*BBCORE*/
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 "bbsavestate", _hoc_bbsavestate,
 "nernst", _hoc_nernst,
 "setRNG", _hoc_setRNG,
 "toggleVerbose", _hoc_toggleVerbose,
 "toggleLTPlasticity", _hoc_toggleLTPlasticity,
 "urand", _hoc_urand,
 0, 0
};
 
#endif /*BBCORE*/
#define bbsavestate bbsavestate_GluSynapse
#define nernst nernst_GluSynapse
#define toggleVerbose toggleVerbose_GluSynapse
#define toggleLTPlasticity toggleLTPlasticity_GluSynapse
#define urand urand_GluSynapse
 inline double bbsavestate( _threadargsproto_ );
 inline double nernst( _threadargsprotocomma_ double , double , double );
 inline double toggleVerbose( _threadargsproto_ );
 inline double toggleLTPlasticity( _threadargsproto_ );
 inline double urand( _threadargsproto_ );
 /* declare global and static user variables */
 
static void _acc_globals_update() {
 }
 
#if 0 /*BBCORE*/
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tau_r_AMPA", "ms",
 "tau_d_AMPA", "ms",
 "E_AMPA", "mV",
 "gmax_AMPA", "nS",
 "tau_r_NMDA", "ms",
 "tau_d_NMDA", "ms",
 "E_NMDA", "mV",
 "Use", "1",
 "Dep", "ms",
 "Fac", "ms",
 "Nrrp", "1",
 "volume_CR", "um3",
 "gca_bar_VDCC", "nS/um2",
 "ljp_VDCC", "mV",
 "vhm_VDCC", "mV",
 "km_VDCC", "mV",
 "vhh_VDCC", "mV",
 "kh_VDCC", "mV",
 "tm_VDCC", "ms",
 "th_VDCC", "ms",
 "gamma_ca_CR", "1",
 "tau_ca_CR", "ms",
 "min_ca_CR", "mM",
 "cao_CR", "mM",
 "tau_GB", "s",
 "theta_d_GB", "mM",
 "theta_p_GB", "mM",
 "gamma_d_GB", "1",
 "gamma_p_GB", "1",
 "rho_star_GB", "1",
 "rho0_GB", "1",
 "enable_GB", "1",
 "tau_Use_GB", "s",
 "Use_d_GB", "1",
 "Use_p_GB", "1",
 "NMDA_ratio", "1",
 "mg", "mM",
 "A_AMPA", "1",
 "B_AMPA", "1",
 "A_NMDA", "1",
 "B_NMDA", "1",
 "m_VDCC", "1",
 "h_VDCC", "1",
 "cai_CR", "mM",
 "Rho_GB", "1",
 "Use_GB", "1",
 "effcai_GB", "1",
 "factor_AMPA", "1",
 "g_AMPA", "uS",
 "i_AMPA", "nA",
 "factor_NMDA", "1",
 "gmax_NMDA", "nS",
 "g_NMDA", "uS",
 "i_NMDA", "nA",
 "u", "1",
 "tsyn", "ms",
 "Psurv", "1",
 "unoccupied", "1",
 "occupied", "1",
 "Pf_NMDA", "1",
 "ica_NMDA", "nA",
 "area_CR", "um2",
 "gca_VDCC", "uS",
 "gca_bar_abs_VDCC", "nS",
 "minf_VDCC", "1",
 "hinf_VDCC", "1",
 "mtau_VDCC", "ms",
 "htau_VDCC", "ms",
 "ica_VDCC", "nA",
 "depress_GB", "1",
 "potentiate_GB", "1",
 "w", "1",
 "g", "uS",
 "i", "nA",
 0,0
};
 
#endif /*BBCORE*/
 static double A_NMDA0 = 0;
 static double A_AMPA0 = 0;
 static double B_NMDA0 = 0;
 static double B_AMPA0 = 0;
 static double Rho_GB0 = 0;
 static double Use_GB0 = 0;
 static double cai_CR0 = 0;
 static double delta_t = 0.01;
 static double effcai_GB0 = 0;
 static double h_VDCC0 = 0;
 static double m_VDCC0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
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
 void _nrn_watch_check(NrnThread*, Memb_list*);
 
#define _watch_array(arg) _ppvar[(arg + 4)*_STRIDE] 
 
#define _nrn_watch_activate(_item)\
  if (_watch_rm == 0) {\
    int _i;\
    for (_i = 1; _i < 5; ++_i) {\
      _watch_array(_i) = 0;\
    }\
    _watch_rm = 1;\
  }\
  _watch_array(_item) = 2 + 
#if 0 /*BBCORE*/
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   Prop* _prop = ((Point_process*)_vptr)->_prop;
   destroy_point_process(_vptr);
}
 
#endif /*BBCORE*/
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "6.2.0",
"GluSynapse",
 "tau_r_AMPA",
 "tau_d_AMPA",
 "E_AMPA",
 "gmax_AMPA",
 "tau_r_NMDA",
 "tau_d_NMDA",
 "E_NMDA",
 "Use",
 "Dep",
 "Fac",
 "Nrrp",
 "volume_CR",
 "gca_bar_VDCC",
 "ljp_VDCC",
 "vhm_VDCC",
 "km_VDCC",
 "vhh_VDCC",
 "kh_VDCC",
 "tm_VDCC",
 "th_VDCC",
 "gamma_ca_CR",
 "tau_ca_CR",
 "min_ca_CR",
 "cao_CR",
 "tau_GB",
 "theta_d_GB",
 "theta_p_GB",
 "gamma_d_GB",
 "gamma_p_GB",
 "rho_star_GB",
 "rho0_GB",
 "enable_GB",
 "tau_Use_GB",
 "Use_d_GB",
 "Use_p_GB",
 "NMDA_ratio",
 "mg",
 "synapseID",
 "verbose",
 "selected_for_report",
 0,
 "factor_AMPA",
 "g_AMPA",
 "i_AMPA",
 "factor_NMDA",
 "gmax_NMDA",
 "g_NMDA",
 "i_NMDA",
 "u",
 "tsyn",
 "Psurv",
 "unoccupied",
 "occupied",
 "Pf_NMDA",
 "ica_NMDA",
 "area_CR",
 "gca_VDCC",
 "gca_bar_abs_VDCC",
 "minf_VDCC",
 "hinf_VDCC",
 "mtau_VDCC",
 "htau_VDCC",
 "ica_VDCC",
 "depress_GB",
 "potentiate_GB",
 "w",
 "g",
 "i",
 0,
 "A_AMPA",
 "B_AMPA",
 "A_NMDA",
 "B_NMDA",
 "m_VDCC",
 "h_VDCC",
 "cai_CR",
 "Rho_GB",
 "Use_GB",
 "effcai_GB",
 0,
 "rng_rel",
 0};
 
static void nrn_alloc(double* _p, Datum* _ppvar, int _type) {
 
#if 0 /*BBCORE*/
 	/*initialize range parameters*/
 	tau_r_AMPA = 0.2;
 	tau_d_AMPA = 1.7;
 	E_AMPA = 0;
 	gmax_AMPA = 1;
 	tau_r_NMDA = 0.29;
 	tau_d_NMDA = 43;
 	E_NMDA = 0;
 	Use = 1;
 	Dep = 100;
 	Fac = 10;
 	Nrrp = 1;
 	volume_CR = 0.087;
 	gca_bar_VDCC = 0.0744;
 	ljp_VDCC = 0;
 	vhm_VDCC = -5.9;
 	km_VDCC = 9.5;
 	vhh_VDCC = -39;
 	kh_VDCC = -9.2;
 	tm_VDCC = 1;
 	th_VDCC = 27;
 	gamma_ca_CR = 0.04;
 	tau_ca_CR = 12;
 	min_ca_CR = 7e-05;
 	cao_CR = 2;
 	tau_GB = 100;
 	theta_d_GB = 0.006;
 	theta_p_GB = 0.001;
 	gamma_d_GB = 100;
 	gamma_p_GB = 450;
 	rho_star_GB = 0.5;
 	rho0_GB = 0;
 	enable_GB = 0;
 	tau_Use_GB = 100;
 	Use_d_GB = 0.2;
 	Use_p_GB = 0.8;
 	NMDA_ratio = 0.71;
 	mg = 1;
 	synapseID = 0;
 	verbose = 0;
 	selected_for_report = 0;
 
#endif /* BBCORE */
 
}
 static void _initlists();
 
#define _tqitem &(_nt->_vdata[_ppvar[3*_STRIDE]])
 
#if NET_RECEIVE_BUFFERING
#undef _tqitem
#define _tqitem _ppvar[3*_STRIDE]
#endif

 static void _net_receive(Point_process*, int, double);
 
#define _psize 91
#define _ppsize 9
 void _GluSynapse_reg() {
	int _vectorized = 1;
  _initlists();
 _mechtype = nrn_get_mechtype(_mechanism[1]);
 if (_mechtype == -1) return;
 _nrn_layout_reg(_mechtype, LAYOUT);
 
#if 0 /*BBCORE*/
 
#endif /*BBCORE*/
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, NULL, nrn_state, nrn_init,
	 hoc_nrnpointerindex,
	 NULL/*_hoc_create_pnt*/, NULL/*_hoc_destroy_pnt*/, /*_member_func,*/
	 1);
   hoc_reg_bbcore_read(_mechtype, bbcore_read);
   hoc_reg_bbcore_write(_mechtype, bbcore_write);
  hoc_register_prop_size(_mechtype, _psize, _ppsize);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "bbcorepointer");
  hoc_register_dparam_semantics(_mechtype, 3, "netsend");
  hoc_register_dparam_semantics(_mechtype, 4, "watch");
  hoc_register_dparam_semantics(_mechtype, 5, "watch");
  hoc_register_dparam_semantics(_mechtype, 6, "watch");
  hoc_register_dparam_semantics(_mechtype, 7, "watch");
  hoc_register_dparam_semantics(_mechtype, 8, "watch");
 hoc_register_watch_check(_nrn_watch_check, _mechtype);
 
#if NET_RECEIVE_BUFFERING
  hoc_register_net_receive_buffering(_net_buf_receive, _mechtype);
#endif
 
#if NET_RECEIVE_BUFFERING
  hoc_register_net_send_buffering(_mechtype);
#endif
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, NULL);
 }
 
double FARADAY = 96485.3;
#pragma acc declare copyin(FARADAY)
 
double PI = 3.14159;
#pragma acc declare copyin(PI)
 
double R = 8.3145;
#pragma acc declare copyin(R)
static char *modelname = "Glutamatergic synapse";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static inline int setRNG(_threadargsproto_);
 
/* _euler_ state _GluSynapse */
#ifndef INSIDE_NMODL
#define INSIDE_NMODL
#endif
 
int _ode_spec1(_threadargsproto_);
/*int _ode_matsol1(_threadargsproto_);*/
 static double *_temp1;
 
#define _slist1 _slist1_GluSynapse
int* _slist1;
#pragma acc declare create(_slist1)

#define _dlist1 _dlist1_GluSynapse
int* _dlist1;
#pragma acc declare create(_dlist1)
 extern int state(_threadargsproto_);
 } using namespace coreneuron; 
/*VERBATIM*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "nrnran123.h"

double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
 namespace coreneuron { 
/*CVODE*/
 int _ode_spec1 (_threadargsproto_) {int _reset = 0; {
   DA_AMPA = - A_AMPA / tau_r_AMPA ;
   DB_AMPA = - B_AMPA / tau_d_AMPA ;
   DA_NMDA = - A_NMDA / tau_r_NMDA ;
   DB_NMDA = - B_NMDA / tau_d_NMDA ;
   Dm_VDCC = ( minf_VDCC - m_VDCC ) / mtau_VDCC ;
   Dh_VDCC = ( hinf_VDCC - h_VDCC ) / htau_VDCC ;
   Dcai_CR = - ( 1e-9 ) * ( ica_NMDA + ica_VDCC ) * gamma_ca_CR / ( ( 1e-15 ) * volume_CR * 2.0 * FARADAY ) - ( cai_CR - min_ca_CR ) / tau_ca_CR ;
   Deffcai_GB = - 0.005 * effcai_GB + ( cai_CR - min_ca_CR ) ;
   DRho_GB = ( - Rho_GB * ( 1.0 - Rho_GB ) * ( rho_star_GB - Rho_GB ) + potentiate_GB * gamma_p_GB * ( 1.0 - Rho_GB ) - depress_GB * gamma_d_GB * Rho_GB ) / ( ( 1e3 ) * tau_GB ) ;
   DUse_GB = ( Use_d_GB + Rho_GB * ( Use_p_GB - Use_d_GB ) - Use_GB ) / ( ( 1e3 ) * tau_Use_GB ) ;
   }
 return _reset;
}
 int _ode_matsol1 (_threadargsproto_) {
 DA_AMPA = DA_AMPA  / (1. - dt*( ( - 1.0 ) / tau_r_AMPA )) ;
 DB_AMPA = DB_AMPA  / (1. - dt*( ( - 1.0 ) / tau_d_AMPA )) ;
 DA_NMDA = DA_NMDA  / (1. - dt*( ( - 1.0 ) / tau_r_NMDA )) ;
 DB_NMDA = DB_NMDA  / (1. - dt*( ( - 1.0 ) / tau_d_NMDA )) ;
 Dm_VDCC = Dm_VDCC  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau_VDCC )) ;
 Dh_VDCC = Dh_VDCC  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau_VDCC )) ;
 Dcai_CR = Dcai_CR  / (1. - dt*( ( - ( ( 1.0 ) ) / tau_ca_CR ) )) ;
 Deffcai_GB = Deffcai_GB  / (1. - dt*( (- 0.005)*(1.0) )) ;
 DRho_GB = DRho_GB  / (1. - dt*( ( ( ((((- 1.0)*(( 1.0 - Rho_GB )) + (- Rho_GB)*(( ( - 1.0 ) ))))*(( rho_star_GB - Rho_GB )) + (- Rho_GB * ( 1.0 - Rho_GB ))*(( ( - 1.0 ) ))) + (potentiate_GB * gamma_p_GB)*(( ( - 1.0 ) )) - (depress_GB * gamma_d_GB)*(1.0) ) ) / ( ( 1e3 ) * tau_GB ) )) ;
 DUse_GB = DUse_GB  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ( ( 1e3 ) * tau_Use_GB ) )) ;
 return 0;
}
 /*END CVODE*/
 
int state (_threadargsproto_) {int _reset=0; int error = 0;
 {
   DA_AMPA = - A_AMPA / tau_r_AMPA ;
   DB_AMPA = - B_AMPA / tau_d_AMPA ;
   DA_NMDA = - A_NMDA / tau_r_NMDA ;
   DB_NMDA = - B_NMDA / tau_d_NMDA ;
   Dm_VDCC = ( minf_VDCC - m_VDCC ) / mtau_VDCC ;
   Dh_VDCC = ( hinf_VDCC - h_VDCC ) / htau_VDCC ;
   Dcai_CR = - ( 1e-9 ) * ( ica_NMDA + ica_VDCC ) * gamma_ca_CR / ( ( 1e-15 ) * volume_CR * 2.0 * FARADAY ) - ( cai_CR - min_ca_CR ) / tau_ca_CR ;
   Deffcai_GB = - 0.005 * effcai_GB + ( cai_CR - min_ca_CR ) ;
   DRho_GB = ( - Rho_GB * ( 1.0 - Rho_GB ) * ( rho_star_GB - Rho_GB ) + potentiate_GB * gamma_p_GB * ( 1.0 - Rho_GB ) - depress_GB * gamma_d_GB * Rho_GB ) / ( ( 1e3 ) * tau_GB ) ;
   DUse_GB = ( Use_d_GB + Rho_GB * ( Use_p_GB - Use_d_GB ) - Use_GB ) / ( ( 1e3 ) * tau_Use_GB ) ;
   }
 return _reset;}
 
#if NET_RECEIVE_BUFFERING 
#undef t
#define t _nrb_t
static inline void _net_receive_kernel(double, Point_process*, int _weight_index, double _flag);
static void _net_buf_receive(NrnThread* _nt) {
  if (!_nt->_ml_list) { return; }
  Memb_list* _ml = _nt->_ml_list[_mechtype];
  if (!_ml) { return; }
  NetReceiveBuffer_t* _nrb = _ml->_net_receive_buffer; 
  int _di;
  int stream_id = _nt->stream_id;
  Point_process* _pnt = _nt->pntprocs;
  int _pnt_length = _nt->n_pntproc - _nrb->_pnt_offset;
  int _displ_cnt = _nrb->_displ_cnt;
  _PRAGMA_FOR_NETRECV_ACC_LOOP_ 
  for (_di = 0; _di < _displ_cnt; ++_di) {
    int _inrb;
    int _di0 = _nrb->_displ[_di];
    int _di1 = _nrb->_displ[_di + 1];
    for (_inrb = _di0; _inrb < _di1; ++_inrb) {
      int _i = _nrb->_nrb_index[_inrb];
      int _j = _nrb->_pnt_index[_i];
      int _k = _nrb->_weight_index[_i];
      double _nrt = _nrb->_nrb_t[_i];
      double _nrflag = _nrb->_nrb_flag[_i];
      _net_receive_kernel(_nrt, _pnt + _j, _k, _nrflag);
    }
  }
  #pragma acc wait(stream_id)
  _nrb->_displ_cnt = 0;
  _nrb->_cnt = 0;
  /*printf("_net_buf_receive__GluSynapse  %d\n", _nt->_id);*/
 
  {
  NetSendBuffer_t* _nsb = _ml->_net_send_buffer;
#if defined(_OPENACC) && !defined(DISABLE_OPENACC)
  #pragma acc update self(_nsb->_cnt) if(_nt->compute_gpu)
  update_net_send_buffer_on_host(_nt, _nsb);
#endif
  int _i;
  for (_i=0; _i < _nsb->_cnt; ++_i) {
    net_sem_from_gpu(_nsb->_sendtype[_i], _nsb->_vdata_index[_i],
      _nsb->_weight_index[_i], _nt->_id, _nsb->_pnt_index[_i],
      _nsb->_nsb_t[_i], _nsb->_nsb_flag[_i]);
  }
  _nsb->_cnt = 0;
#if defined(_OPENACC) && !defined(DISABLE_OPENACC)
  #pragma acc update device(_nsb->_cnt) if (_nt->compute_gpu)
#endif
  }
 
}
 
static void _net_send_buffering(NetSendBuffer_t* _nsb, int _sendtype, int _i_vdata, int _weight_index,
 int _ipnt, double _t, double _flag) {
  int _i = 0;
  #pragma acc atomic capture
  _i = _nsb->_cnt++;
  if (_i >= _nsb->_size) {
  }
  _nsb->_sendtype[_i] = _sendtype;
  _nsb->_vdata_index[_i] = _i_vdata;
  _nsb->_weight_index[_i] = _weight_index;
  _nsb->_pnt_index[_i] = _ipnt;
  _nsb->_nsb_t[_i] = _t;
  _nsb->_nsb_flag[_i] = _flag;
}
 
static void _net_receive (Point_process* _pnt, int _weight_index, double _lflag) {
  NrnThread* _nt = nrn_threads + _pnt->_tid;
  NetReceiveBuffer_t* _nrb = _nt->_ml_list[_mechtype]->_net_receive_buffer;
  if (_nrb->_cnt >= _nrb->_size){
    realloc_net_receive_buffer(_nt, _nt->_ml_list[_mechtype]);
  }
  _nrb->_pnt_index[_nrb->_cnt] = _pnt - _nt->pntprocs;
  _nrb->_weight_index[_nrb->_cnt] = _weight_index;
  _nrb->_nrb_t[_nrb->_cnt] = _nt->_t;
  _nrb->_nrb_flag[_nrb->_cnt] = _lflag;
  ++_nrb->_cnt;
}
 
static void _net_receive_kernel(double _nrb_t, Point_process* _pnt, int _weight_index, double _lflag)
#else
 
static void _net_receive (Point_process* _pnt, int _weight_index, double _lflag) 
#endif
 
{  double* _p; Datum* _ppvar; ThreadDatum* _thread; double v;
   Memb_list* _ml; int _cntml_padded, _cntml_actual; int _iml; double* _args;
   int _watch_rm = 0;
 
   NrnThread* _nt;
   int _tid = _pnt->_tid; 
   _nt = nrn_threads + _tid;
   _thread = (ThreadDatum*)0; 
   double *_weights = _nt->_weights;
   _args = _weights + _weight_index;
   _ml = _nt->_ml_list[_pnt->_type];
   _cntml_actual = _ml->_nodecount;
   _cntml_padded = _ml->_nodecount_padded;
   _iml = _pnt->_i_instance;
#if LAYOUT == 1 /*AoS*/
   _p = _ml->_data + _iml*_psize; _ppvar = _ml->_pdata + _iml*_ppsize;
#endif
#if LAYOUT == 0 /*SoA*/
   _p = _ml->_data; _ppvar = _ml->_pdata;
#endif
#if LAYOUT > 1 /*AoSoA*/
#error AoSoA not implemented.
#endif
  #if !defined(_OPENACC) 
 assert(_tsav <= t); 
 #endif 
 _tsav = t; 
#if !NET_RECEIVE_BUFFERING
  if (_lflag == 1. ) {*(_tqitem) = 0;}
#endif
 {
   double _lresult , _lves , _loccu , _lUse_actual ;
 if ( _lflag  == 1.0 ) {
     _nrn_watch_activate(1)  ( effcai_GB > theta_d_GB ) ; /* 2.0 */
 _nrn_watch_activate(2)  ( effcai_GB < theta_d_GB ) ; /* 3.0 */
 _nrn_watch_activate(3)  ( effcai_GB > theta_p_GB ) ; /* 4.0 */
 _nrn_watch_activate(4)  ( effcai_GB < theta_p_GB ) ; /* 5.0 */
 if ( verbose > 0.0 ) {
        printf ( "Flag 1, Initialize watch calls\n" ) ;
        }
     }
   else if ( _lflag  == 2.0 ) {
     depress_GB = 1.0 ;
     }
   else if ( _lflag  == 3.0 ) {
     depress_GB = 0.0 ;
     }
   else if ( _lflag  == 4.0 ) {
     potentiate_GB = 1.0 ;
     }
   else if ( _lflag  == 5.0 ) {
     potentiate_GB = 0.0 ;
     }
   else {
     if (  ! ( _args[0] > 0.0 ) ) {
       
/*VERBATIM*/
            return;
 }
     w = _args[0] ;
     if ( verbose > 0.0 ) {
        printf ( "t = %g, incoming spike at synapse %g\n" , t , synapseID ) ;
        }
     if ( enable_GB  == 1.0 ) {
       _lUse_actual = Use_GB ;
       }
     else {
       _lUse_actual = Use ;
       }
     if ( Fac > 0.0 ) {
       u = u * exp ( - ( t - tsyn ) / Fac ) ;
       u = u + _lUse_actual * ( 1.0 - u ) ;
       }
     else {
       u = _lUse_actual ;
       }
     {int  _lcounter ;for ( _lcounter = 0 ; _lcounter <= ( ((int) unoccupied ) - 1 ) ; _lcounter ++ ) {
       Psurv = exp ( - ( t - tsyn ) / Dep ) ;
       _lresult = urand ( _threadargs_ ) ;
       if ( _lresult > Psurv ) {
         occupied = occupied + 1.0 ;
         if ( verbose > 0.0 ) {
            printf ( "\tRecovered 1 vesicle, P = %g R = %g\n" , Psurv , _lresult ) ;
            }
         }
       } }
     _lves = 0.0 ;
     _loccu = occupied - 1.0 ;
     {int  _lcounter ;for ( _lcounter = 0 ; _lcounter <= ((int) _loccu ) ; _lcounter ++ ) {
       _lresult = urand ( _threadargs_ ) ;
       if ( _lresult < u ) {
         occupied = occupied - 1.0 ;
         _lves = _lves + 1.0 ;
         }
       } }
     unoccupied = Nrrp - occupied ;
     tsyn = t ;
     A_AMPA = A_AMPA + _lves / Nrrp * factor_AMPA ;
     B_AMPA = B_AMPA + _lves / Nrrp * factor_AMPA ;
     A_NMDA = A_NMDA + _lves / Nrrp * factor_NMDA ;
     B_NMDA = B_NMDA + _lves / Nrrp * factor_NMDA ;
     if ( verbose > 0.0 ) {
        printf ( "\tReleased %g vesicles out of %g\n" , _lves , Nrrp ) ;
        }
     }
   } 
#if NET_RECEIVE_BUFFERING
#undef t
#define t _nt->_t
#endif
 }
 
double nernst ( _threadargsprotocomma_ double _lci , double _lco , double _lz ) {
   double _lnernst;
 _lnernst = ( 1000.0 ) * R * ( celsius + 273.15 ) / ( _lz * FARADAY ) * log ( _lco / _lci ) ;
   if ( verbose > 1.0 ) {
      printf ( "nernst:%f R:%f celsius:%f \n" , _lnernst , R , celsius ) ;
      }
   
return _lnernst;
 }
 
#if 0 /*BBCORE*/
 
static double _hoc_nernst(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  nernst ( _threadargs_, *getarg(1) , *getarg(2) , *getarg(3) );
 return(_r);
}
 
#endif /*BBCORE*/
 
static int  setRNG ( _threadargsproto_ ) {
   
/*VERBATIM*/
    #ifndef CORENEURON_BUILD
    // For compatibility, allow for either MCellRan4 or Random123
    // Distinguish by the arg types
    // Object => MCellRan4, seeds (double) => Random123
    usingR123 = 0;
    if( ifarg(1) && hoc_is_double_arg(1) ) {
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng_rel);
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
    } else if( ifarg(1) ) {   // not a double, so assume hoc object type
        void** pv = (void**)(&_p_rng_rel);
        *pv = nrn_random_arg(1);
    } else {  // no arg, so clear pointer
        void** pv = (void**)(&_p_rng_rel);
        *pv = (void*)0;
    }
    #endif
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_setRNG(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 setRNG ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
double urand ( _threadargsproto_ ) {
   double _lurand;
 
/*VERBATIM*/
    double value;
    if ( usingR123 ) {
        value = nrnran123_dblpick((nrnran123_State*)_p_rng_rel);
    } else if (_p_rng_rel) {
        #ifndef CORENEURON_BUILD
        value = nrn_random_pick(_p_rng_rel);
        #endif
    } else {
        value = 0.0;
    }
    _lurand = value;
 
return _lurand;
 }
 
#if 0 /*BBCORE*/
 
static double _hoc_urand(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  urand ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
double toggleLTPlasticity ( _threadargsproto_ ) {
   double _ltoggleLTPlasticity;
 enable_GB = 1.0 - enable_GB ;
   
return _ltoggleLTPlasticity;
 }
 
#if 0 /*BBCORE*/
 
static double _hoc_toggleLTPlasticity(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  toggleLTPlasticity ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
double toggleVerbose ( _threadargsproto_ ) {
   double _ltoggleVerbose;
 verbose = 1.0 - verbose ;
   
return _ltoggleVerbose;
 }
 
#if 0 /*BBCORE*/
 
static double _hoc_toggleVerbose(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  toggleVerbose ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
double bbsavestate ( _threadargsproto_ ) {
   double _lbbsavestate;
 _lbbsavestate = 0.0 ;
   
/*VERBATIM*/
#ifndef CORENEURON_BUILD
        /* first arg is direction (0 save, 1 restore), second is array*/
        /* if first arg is -1, fill xdir with the size of the array */
        double *xdir, *xval, *hoc_pgetarg();
        long nrn_get_random_sequence(void* r);
        void nrn_set_random_sequence(void* r, int val);
        xdir = hoc_pgetarg(1);
        xval = hoc_pgetarg(2);
        if (_p_rng_rel) {
            // tell how many items need saving
            if (*xdir == -1) {  // count items
                if( usingR123 ) {
                    *xdir = 2.0;
                } else {
                    *xdir = 1.0;
                }
                return 0.0;
            } else if(*xdir ==0 ) {  // save
                if( usingR123 ) {
                    uint32_t seq;
                    unsigned char which;
                    nrnran123_getseq( (nrnran123_State*)_p_rng_rel, &seq, &which );
                    xval[0] = (double) seq;
                    xval[1] = (double) which;
                } else {
                    xval[0] = (double)nrn_get_random_sequence(_p_rng_rel);
                }
            } else {  // restore
                if( usingR123 ) {
                    nrnran123_setseq( (nrnran123_State*)_p_rng_rel, (uint32_t)xval[0], (char)xval[1] );
                } else {
                    nrn_set_random_sequence(_p_rng_rel, (long)(xval[0]));
                }
            }
        }
#endif
 
return _lbbsavestate;
 }
 
#if 0 /*BBCORE*/
 
static double _hoc_bbsavestate(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  bbsavestate ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 } using namespace coreneuron; 
/*VERBATIM*/

static void bbcore_write(double* dArray, int* iArray, int* doffset, int* ioffset, _threadargsproto_) {
    // make sure offset array non-null
    if (iArray) {

        // get handle to random123 instance
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng_rel);

        // get location for storing ids
        uint32_t* ia = ((uint32_t*)iArray) + *ioffset;

        // retrieve/store identifier seeds
        nrnran123_getids3(*pv, ia, ia+1, ia+2);

        // retrieve/store stream sequence
        unsigned char which;
        nrnran123_getseq(*pv, ia+3, &which);
        ia[4] = (int)which;
    }

    // increment integer offset (2 identifier), no double data
    *ioffset += 5;
    *doffset += 0;
}

static void bbcore_read(double* dArray, int* iArray, int* doffset, int* ioffset, _threadargsproto_) {
    // make sure it's not previously set
    assert(!_p_rng_rel);

    uint32_t* ia = ((uint32_t*)iArray) + *ioffset;

    // make sure non-zero identifier seeds
    if (ia[0] != 0 || ia[1] != 0 || ia[2] != 0) {
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng_rel);

        // get new stream
        *pv = nrnran123_newstream3(ia[0], ia[1], ia[2]);

        // restore sequence
        nrnran123_setseq(*pv, ia[3], (char)ia[4]);
    }

    // increment intger offset (2 identifiers), no double data
    *ioffset += 5;
    *doffset += 0;
}

 namespace coreneuron { 
void _nrn_watch_check(NrnThread* _nt, Memb_list* _ml) {
  double* _p; Datum* _ppvar; ThreadDatum* _thread;
  int* _ni; double v; int _iml, _cntml_padded, _cntml_actual;
  _cntml_actual = _ml->_nodecount;
  _cntml_padded = _ml->_nodecount_padded;
  _ni = _ml->_nodeindices;
  _thread = _ml->_thread;
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
    v = _vec_v[_ni[_iml]];
 
    if (_watch_array(1)&2) {
      if  ( effcai_GB > theta_d_GB ) {
        if ((_watch_array(1)&1) == 0) {
          #if NET_RECEIVE_BUFFERING
          _net_send_buffering(_ml->_net_send_buffer, 0, _tqitem, 0, _ppvar[1*_STRIDE], t +  0.0 , 2.0 );
          #else
          net_send ( _tqitem, -1, (Point_process*) _nt->_vdata[_ppvar[1*_STRIDE]], t +  0.0 , 2.0 ) ;
          #endif
        }
        _watch_array(1) = 3;
      }else{
        _watch_array(1) = 2;
      }
    }
 
    if (_watch_array(2)&2) {
      if  ( effcai_GB < theta_d_GB ) {
        if ((_watch_array(2)&1) == 0) {
          #if NET_RECEIVE_BUFFERING
          _net_send_buffering(_ml->_net_send_buffer, 0, _tqitem, 0, _ppvar[1*_STRIDE], t +  0.0 , 3.0 );
          #else
          net_send ( _tqitem, -1, (Point_process*) _nt->_vdata[_ppvar[1*_STRIDE]], t +  0.0 , 3.0 ) ;
          #endif
        }
        _watch_array(2) = 3;
      }else{
        _watch_array(2) = 2;
      }
    }
 
    if (_watch_array(3)&2) {
      if  ( effcai_GB > theta_p_GB ) {
        if ((_watch_array(3)&1) == 0) {
          #if NET_RECEIVE_BUFFERING
          _net_send_buffering(_ml->_net_send_buffer, 0, _tqitem, 0, _ppvar[1*_STRIDE], t +  0.0 , 4.0 );
          #else
          net_send ( _tqitem, -1, (Point_process*) _nt->_vdata[_ppvar[1*_STRIDE]], t +  0.0 , 4.0 ) ;
          #endif
        }
        _watch_array(3) = 3;
      }else{
        _watch_array(3) = 2;
      }
    }
 
    if (_watch_array(4)&2) {
      if  ( effcai_GB < theta_p_GB ) {
        if ((_watch_array(4)&1) == 0) {
          #if NET_RECEIVE_BUFFERING
          _net_send_buffering(_ml->_net_send_buffer, 0, _tqitem, 0, _ppvar[1*_STRIDE], t +  0.0 , 5.0 );
          #else
          net_send ( _tqitem, -1, (Point_process*) _nt->_vdata[_ppvar[1*_STRIDE]], t +  0.0 , 5.0 ) ;
          #endif
        }
        _watch_array(4) = 3;
      }else{
        _watch_array(4) = 2;
      }
    }
   }
 
#if NET_RECEIVE_BUFFERING
  NetSendBuffer_t* _nsb = _ml->_net_send_buffer;
#if defined(_OPENACC) && !defined(DISABLE_OPENACC)
  #pragma acc wait(stream_id)
  #pragma acc update self(_nsb->_cnt) if(_nt->compute_gpu)
  update_net_send_buffer_on_host(_nt, _nsb);
#endif
  {int _i;
  for (_i=0; _i < _nsb->_cnt; ++_i) {
    net_sem_from_gpu(_nsb->_sendtype[_i], _nsb->_vdata_index[_i],
      _nsb->_weight_index[_i], _nt->_id, _nsb->_pnt_index[_i],
      _nsb->_nsb_t[_i], _nsb->_nsb_flag[_i]);
  }}
  _nsb->_cnt = 0;
#if defined(_OPENACC) && !defined(DISABLE_OPENACC)
  #pragma acc update device(_nsb->_cnt) if(_nt->compute_gpu)
#endif
#endif
 }

static inline void initmodel(_threadargsproto_) {
  int _i; double _save;{
  Memb_list* _ml = _nt->_ml_list[_mechtype];
  A_NMDA = A_NMDA0;
  A_AMPA = A_AMPA0;
  B_NMDA = B_NMDA0;
  B_AMPA = B_AMPA0;
  Rho_GB = Rho_GB0;
  Use_GB = Use_GB0;
  cai_CR = cai_CR0;
  effcai_GB = effcai_GB0;
  h_VDCC = h_VDCC0;
  m_VDCC = m_VDCC0;
 {
   double _ltp_AMPA , _ltp_NMDA ;
 A_AMPA = 0.0 ;
   B_AMPA = 0.0 ;
   _ltp_AMPA = ( tau_r_AMPA * tau_d_AMPA ) / ( tau_d_AMPA - tau_r_AMPA ) * log ( tau_d_AMPA / tau_r_AMPA ) ;
   factor_AMPA = 1.0 / ( - exp ( - _ltp_AMPA / tau_r_AMPA ) + exp ( - _ltp_AMPA / tau_d_AMPA ) ) ;
   A_NMDA = 0.0 ;
   B_NMDA = 0.0 ;
   _ltp_NMDA = ( tau_r_NMDA * tau_d_NMDA ) / ( tau_d_NMDA - tau_r_NMDA ) * log ( tau_d_NMDA / tau_r_NMDA ) ;
   factor_NMDA = 1.0 / ( - exp ( - _ltp_NMDA / tau_r_NMDA ) + exp ( - _ltp_NMDA / tau_d_NMDA ) ) ;
   gmax_NMDA = gmax_AMPA * NMDA_ratio ;
   tsyn = 0.0 ;
   Psurv = 0.0 ;
   u = 0.0 ;
   unoccupied = 0.0 ;
   occupied = Nrrp ;
   Pf_NMDA = ( 4.0 * cao_CR ) / ( 4.0 * cao_CR + ( 1.0 / 1.38 ) * 120.0 ) * 0.6 ;
    area_CR = 4.0 * PI * pow( ( 3.0 / 4.0 * volume_CR * 1.0 / PI ) , ( 2.0 / 3.0 ) ) ;
    gca_bar_abs_VDCC = gca_bar_VDCC * area_CR ;
   mtau_VDCC = tm_VDCC ;
   htau_VDCC = th_VDCC ;
   cai_CR = min_ca_CR ;
   Rho_GB = rho0_GB ;
   Use_GB = Use ;
   effcai_GB = 0.0 ;
   depress_GB = 0.0 ;
   potentiate_GB = 0.0 ;
   w = 0.0 ;
   
#if NET_RECEIVE_BUFFERING
    _net_send_buffering(_ml->_net_send_buffer, 0, _tqitem, 0,  _ppvar[1*_STRIDE], t +  0.0 , 1.0 );
#else
 net_send ( _tqitem, -1, (Point_process*) _nt->_vdata[_ppvar[1*_STRIDE]], t +  0.0 , 1.0 ) ;
   
#endif
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
  #pragma acc update device (_mechtype) if(_nt->compute_gpu)

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
 _tsav = -1e20;
    _v = _vec_v[_nd_idx];
    _PRCELLSTATE_V
 v = _v;
 _PRCELLSTATE_V
 initmodel(_threadargs_);
}
  }

#if NET_RECEIVE_BUFFERING
  NetSendBuffer_t* _nsb = _ml->_net_send_buffer;
#if defined(_OPENACC) && !defined(DISABLE_OPENACC)
  #pragma acc wait(stream_id)
  #pragma acc update self(_nsb->_cnt) if(_nt->compute_gpu)
  update_net_send_buffer_on_host(_nt, _nsb);
#endif
  {int _i;
  for (_i=0; _i < _nsb->_cnt; ++_i) {
    net_sem_from_gpu(_nsb->_sendtype[_i], _nsb->_vdata_index[_i],
      _nsb->_weight_index[_i], _nt->_id, _nsb->_pnt_index[_i],
      _nsb->_nsb_t[_i], _nsb->_nsb_flag[_i]);
  }}
  _nsb->_cnt = 0;
#if defined(_OPENACC) && !defined(DISABLE_OPENACC)
  #pragma acc update device(_nsb->_cnt) if(_nt->compute_gpu)
#endif
#endif
}

static double _nrn_current(_threadargsproto_, double _v){double _current=0.;v=_v;{ {
   double _lEca_syn , _lmggate ;
 g_AMPA = w * ( 1e-3 ) * gmax_AMPA * ( B_AMPA - A_AMPA ) ;
   i_AMPA = g_AMPA * ( v - E_AMPA ) ;
   _lmggate = 1.0 / ( 1.0 + exp ( 0.062 * - ( v ) ) * ( mg / 3.57 ) ) ;
   g_NMDA = w * ( 1e-3 ) * gmax_NMDA * _lmggate * ( B_NMDA - A_NMDA ) ;
   i_NMDA = g_NMDA * ( v - E_NMDA ) ;
   ica_NMDA = Pf_NMDA * g_NMDA * ( v - 40.0 ) ;
   _lEca_syn = nernst ( _threadargscomma_ cai_CR , cao_CR , 2.0 ) ;
   gca_VDCC = ( 1e-3 ) * gca_bar_abs_VDCC * m_VDCC * m_VDCC * h_VDCC ;
   ica_VDCC = gca_VDCC * ( v - _lEca_syn ) ;
   minf_VDCC = 1.0 / ( 1.0 + exp ( ( ( vhm_VDCC - ljp_VDCC ) - v ) / km_VDCC ) ) ;
   hinf_VDCC = 1.0 / ( 1.0 + exp ( ( ( vhh_VDCC - ljp_VDCC ) - v ) / kh_VDCC ) ) ;
   g = g_AMPA + g_NMDA ;
   i = i_AMPA + i_NMDA + ica_VDCC ;
   }
 _current += i;

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
double * _vec_shadow_rhs = _nt->_shadow_rhs;
double * _vec_shadow_d = _nt->_shadow_d;

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
_PRAGMA_FOR_CUR_SYN_ACC_LOOP_
for (_iml = 0; _iml < _cntml_actual; ++_iml) {
#else /* LAYOUT > 1 */ /*AoSoA*/
#error AoSoA not implemented.
for (;;) { /* help clang-format properly indent */
#endif
    int _nd_idx = _ni[_iml];
    _v = _vec_v[_nd_idx];
    _PRCELLSTATE_V
 _g = _nrn_current(_threadargs_, _v + .001);
 	{ _rhs = _nrn_current(_threadargs_, _v);
 	}
 _g = (_g - _rhs)/.001;
 double _mfact =  1.e2/(_nd_area);
 _g *=  _mfact;
 _rhs *= _mfact;
 _PRCELLSTATE_G


#ifdef _OPENACC
  if(_nt->compute_gpu) {
    #pragma acc atomic update
    _vec_rhs[_nd_idx] -= _rhs;
    #pragma acc atomic update
    _vec_d[_nd_idx] += _g;
  } else {
    _vec_shadow_rhs[_iml] = _rhs;
    _vec_shadow_d[_iml] = _g;
  }
#else
  _vec_shadow_rhs[_iml] = _rhs;
  _vec_shadow_d[_iml] = _g;
#endif
 }
#ifdef _OPENACC
    if(!(_nt->compute_gpu)) { 
        for (_iml = 0; _iml < _cntml_actual; ++_iml) {
           int _nd_idx = _ni[_iml];
           _vec_rhs[_nd_idx] -= _vec_shadow_rhs[_iml];
           _vec_d[_nd_idx] += _vec_shadow_d[_iml];
        }
#else
 for (_iml = 0; _iml < _cntml_actual; ++_iml) {
   int _nd_idx = _ni[_iml];
   _vec_rhs[_nd_idx] -= _vec_shadow_rhs[_iml];
   _vec_d[_nd_idx] += _vec_shadow_d[_iml];
#endif
 
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
 {  
  #if !defined(_euler_state_GluSynapse)
    #define _euler_state_GluSynapse 0
  #endif
  euler_thread(10, _slist1, _dlist1, _euler_state_GluSynapse, _threadargs_);
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
 int _cntml_actual=1;
 int _cntml_padded=1;
 int _iml=0;
  if (!_first) return;
 
 _slist1 = (int*)malloc(sizeof(int)*10);
 _dlist1 = (int*)malloc(sizeof(int)*10);
 _slist1[0] = &(A_AMPA) - _p;  _dlist1[0] = &(DA_AMPA) - _p;
 _slist1[1] = &(B_AMPA) - _p;  _dlist1[1] = &(DB_AMPA) - _p;
 _slist1[2] = &(A_NMDA) - _p;  _dlist1[2] = &(DA_NMDA) - _p;
 _slist1[3] = &(B_NMDA) - _p;  _dlist1[3] = &(DB_NMDA) - _p;
 _slist1[4] = &(m_VDCC) - _p;  _dlist1[4] = &(Dm_VDCC) - _p;
 _slist1[5] = &(h_VDCC) - _p;  _dlist1[5] = &(Dh_VDCC) - _p;
 _slist1[6] = &(cai_CR) - _p;  _dlist1[6] = &(Dcai_CR) - _p;
 _slist1[7] = &(effcai_GB) - _p;  _dlist1[7] = &(Deffcai_GB) - _p;
 _slist1[8] = &(Rho_GB) - _p;  _dlist1[8] = &(DRho_GB) - _p;
 _slist1[9] = &(Use_GB) - _p;  _dlist1[9] = &(DUse_GB) - _p;
 #pragma acc enter data copyin(_slist1[0:10])
 #pragma acc enter data copyin(_dlist1[0:10])

_first = 0;
}
} // namespace coreneuron_lib
