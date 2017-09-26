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
 
#define _thread_present_ /**/ , _slist1[0:4], _dlist1[0:4] 
 
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
 
#define nrn_init _nrn_init__ProbAMPANMDA_EMS
#define nrn_cur _nrn_cur__ProbAMPANMDA_EMS
#define nrn_cur_parallel _nrn_cur_parallel__ProbAMPANMDA_EMS
#define _nrn_current _nrn_current__ProbAMPANMDA_EMS
#define nrn_jacob _nrn_jacob__ProbAMPANMDA_EMS
#define nrn_state _nrn_state__ProbAMPANMDA_EMS
#define initmodel initmodel__ProbAMPANMDA_EMS
#define _net_receive _net_receive__ProbAMPANMDA_EMS
#define _net_receive2 _net_receive2__ProbAMPANMDA_EMS
#define nrn_state_launcher nrn_state_ProbAMPANMDA_EMS_launcher
#define nrn_cur_launcher nrn_cur_ProbAMPANMDA_EMS_launcher
#define nrn_jacob_launcher nrn_jacob_ProbAMPANMDA_EMS_launcher 
#if NET_RECEIVE_BUFFERING
#define _net_buf_receive _net_buf_receive_ProbAMPANMDA_EMS
static void _net_buf_receive(_NrnThread*);
#endif
 
#define setRNG setRNG_ProbAMPANMDA_EMS 
#define state state_ProbAMPANMDA_EMS 
#define _ode_matsol1 _nrn_ode_matsol1__ProbAMPANMDA_EMS
#define _ode_spec1 _nrn_ode_spec1__ProbAMPANMDA_EMS

 
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
#define tau_r_AMPA _p[0*_STRIDE]
#define tau_d_AMPA _p[1*_STRIDE]
#define tau_r_NMDA _p[2*_STRIDE]
#define tau_d_NMDA _p[3*_STRIDE]
#define Use _p[4*_STRIDE]
#define Dep _p[5*_STRIDE]
#define Fac _p[6*_STRIDE]
#define e _p[7*_STRIDE]
#define mg _p[8*_STRIDE]
#define u0 _p[9*_STRIDE]
#define synapseID _p[10*_STRIDE]
#define verboseLevel _p[11*_STRIDE]
#define NMDA_ratio _p[12*_STRIDE]
#define i _p[13*_STRIDE]
#define i_AMPA _p[14*_STRIDE]
#define i_NMDA _p[15*_STRIDE]
#define g_AMPA _p[16*_STRIDE]
#define g_NMDA _p[17*_STRIDE]
#define g _p[18*_STRIDE]
#define A_AMPA_step _p[19*_STRIDE]
#define B_AMPA_step _p[20*_STRIDE]
#define A_NMDA_step _p[21*_STRIDE]
#define B_NMDA_step _p[22*_STRIDE]
#define Rstate _p[23*_STRIDE]
#define tsyn_fac _p[24*_STRIDE]
#define tsyn _p[25*_STRIDE]
#define u _p[26*_STRIDE]
#define A_AMPA _p[27*_STRIDE]
#define B_AMPA _p[28*_STRIDE]
#define A_NMDA _p[29*_STRIDE]
#define B_NMDA _p[30*_STRIDE]
#define factor_AMPA _p[31*_STRIDE]
#define factor_NMDA _p[32*_STRIDE]
#define mggate _p[33*_STRIDE]
#define DA_AMPA _p[34*_STRIDE]
#define DB_AMPA _p[35*_STRIDE]
#define DA_NMDA _p[36*_STRIDE]
#define DB_NMDA _p[37*_STRIDE]
#define _v_unused _p[38*_STRIDE]
#define _g_unused _p[39*_STRIDE]
#define _tsav _p[40*_STRIDE]
 
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
#define _p_rng	_nt->_vdata[_ppvar[2*_STRIDE]]
 
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
 static int hoc_nrnpointerindex =  2;
 static ThreadDatum* _extcall_thread;
 /* external NEURON variables */
 
#if 0 /*BBCORE*/
 /* declaration of user functions */
 static double _hoc_setRNG();
 static double _hoc_toggleVerbose();
 static double _hoc_urand();
 
#endif /*BBCORE*/
 static int _mechtype;
 extern int nrn_get_mechtype();
extern void hoc_register_prop_size(int, int, int);
extern Memb_func* memb_func;
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
 "setRNG", _hoc_setRNG,
 "toggleVerbose", _hoc_toggleVerbose,
 "urand", _hoc_urand,
 0, 0
};
 
#endif /*BBCORE*/
#define toggleVerbose toggleVerbose_ProbAMPANMDA_EMS
#define urand urand_ProbAMPANMDA_EMS
 inline double toggleVerbose( _threadargsproto_ );
 inline double urand( _threadargsproto_ );
 /* declare global and static user variables */
#define gmax gmax_ProbAMPANMDA_EMS
 double gmax = 0.001;
 #pragma acc declare copyin (gmax)
 
static void _acc_globals_update() {
 #pragma acc update device (gmax) if(nrn_threads->compute_gpu)
 }
 
#if 0 /*BBCORE*/
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gmax_ProbAMPANMDA_EMS", "uS",
 "tau_r_AMPA", "ms",
 "tau_d_AMPA", "ms",
 "tau_r_NMDA", "ms",
 "tau_d_NMDA", "ms",
 "Use", "1",
 "Dep", "ms",
 "Fac", "ms",
 "e", "mV",
 "mg", "mM",
 "NMDA_ratio", "1",
 "i", "nA",
 "i_AMPA", "nA",
 "i_NMDA", "nA",
 "g_AMPA", "uS",
 "g_NMDA", "uS",
 "g", "uS",
 "Rstate", "1",
 "tsyn_fac", "ms",
 "tsyn", "ms",
 "u", "1",
 0,0
};
 
#endif /*BBCORE*/
 static double A_NMDA0 = 0;
 static double A_AMPA0 = 0;
 static double B_NMDA0 = 0;
 static double B_AMPA0 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "gmax_ProbAMPANMDA_EMS", &gmax_ProbAMPANMDA_EMS,
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
 
#if 0 /*BBCORE*/
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 
#endif /*BBCORE*/
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "6.2.0",
"ProbAMPANMDA_EMS",
 "tau_r_AMPA",
 "tau_d_AMPA",
 "tau_r_NMDA",
 "tau_d_NMDA",
 "Use",
 "Dep",
 "Fac",
 "e",
 "mg",
 "u0",
 "synapseID",
 "verboseLevel",
 "NMDA_ratio",
 0,
 "i",
 "i_AMPA",
 "i_NMDA",
 "g_AMPA",
 "g_NMDA",
 "g",
 "A_AMPA_step",
 "B_AMPA_step",
 "A_NMDA_step",
 "B_NMDA_step",
 "Rstate",
 "tsyn_fac",
 "tsyn",
 "u",
 0,
 "A_AMPA",
 "B_AMPA",
 "A_NMDA",
 "B_NMDA",
 0,
 "rng",
 0};
 

 void _nrn_ode_state_vars__ProbAMPANMDA_EMS(short * count, short** var_offsets, short ** dv_offsets)
 {
     *count = 4;
     *var_offsets = (short*) malloc(sizeof(short)* *count);
     *dv_offsets = (short*) malloc(sizeof(short)* *count);
     *var_offsets[0] = 27;
     *var_offsets[1] = 28;
     *var_offsets[2] = 29;
     *var_offsets[3] = 30;
     *dv_offsets[0] = 34;
     *dv_offsets[1] = 35;
     *dv_offsets[2] = 36;
     *dv_offsets[3] = 37;
 }

static void nrn_alloc(double* _p, Datum* _ppvar, int _type) {
 
#if 0 /*BBCORE*/
 	/*initialize range parameters*/
 	tau_r_AMPA = 0.2;
 	tau_d_AMPA = 1.7;
 	tau_r_NMDA = 0.29;
 	tau_d_NMDA = 43;
 	Use = 1;
 	Dep = 100;
 	Fac = 10;
 	e = 0;
 	mg = 1;
 	u0 = 0;
 	synapseID = 0;
 	verboseLevel = 0;
 	NMDA_ratio = 0.71;
 
#endif /* BBCORE */
 
}
 static void _initlists();
 static void _net_receive(Point_process*, int, double);
 static void _net_init(Point_process*, int, double);
 
#define _psize 41
#define _ppsize 3
 static void bbcore_read(double *, int*, int*, int*, _threadargsproto_);
 extern void hoc_reg_bbcore_read(int, void(*)(double *, int*, int*, int*, _threadargsproto_));
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*f)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(_threadargsproto_, int));
extern void _cvode_abstol( Symbol**, double*, int);

 void _ProbAMPANMDA_EMS_reg() {
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
  hoc_register_prop_size(_mechtype, _psize, _ppsize);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "bbcorepointer");
 
#if NET_RECEIVE_BUFFERING
  hoc_register_net_receive_buffering(_net_buf_receive, _mechtype);
#endif
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_init[_mechtype] = _net_init;
 pnt_receive_size[_mechtype] = 4;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, NULL);
 }
static char *modelname = "Probabilistic AMPA and NMDA receptor with presynaptic short-term plasticity ";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int setRNG(_threadargsproto_);
 
int _ode_spec1(_threadargsproto_);
/*int _ode_matsol1(_threadargsproto_);*/
 
#define _slist1 _slist1_ProbAMPANMDA_EMS
int* _slist1;
#pragma acc declare create(_slist1)

#define _dlist1 _dlist1_ProbAMPANMDA_EMS
int* _dlist1;
#pragma acc declare create(_dlist1)
 static inline int state(_threadargsproto_);
 
/*VERBATIM*/
#include "nrnran123.h"
 
/*CVODE*/
int _ode_spec1 (_threadargsproto_) {int _reset = 0; {
   DA_AMPA = - A_AMPA / tau_r_AMPA ;
   DB_AMPA = - B_AMPA / tau_d_AMPA ;
   DA_NMDA = - A_NMDA / tau_r_NMDA ;
   DB_NMDA = - B_NMDA / tau_d_NMDA ;
   }
 return _reset;
}
int _ode_matsol1 (_threadargsproto_) {
 DA_AMPA = DA_AMPA  / (1. - dt*( ( - 1.0 ) / tau_r_AMPA )) ;
 DB_AMPA = DB_AMPA  / (1. - dt*( ( - 1.0 ) / tau_d_AMPA )) ;
 DA_NMDA = DA_NMDA  / (1. - dt*( ( - 1.0 ) / tau_r_NMDA )) ;
 DB_NMDA = DB_NMDA  / (1. - dt*( ( - 1.0 ) / tau_d_NMDA )) ;
 return 0;
}
 /*END CVODE*/
 static int state (_threadargsproto_) { {
    A_AMPA = A_AMPA + (1. - exp(dt*(( - 1.0 ) / tau_r_AMPA)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau_r_AMPA ) - A_AMPA) ;
    B_AMPA = B_AMPA + (1. - exp(dt*(( - 1.0 ) / tau_d_AMPA)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau_d_AMPA ) - B_AMPA) ;
    A_NMDA = A_NMDA + (1. - exp(dt*(( - 1.0 ) / tau_r_NMDA)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau_r_NMDA ) - A_NMDA) ;
    B_NMDA = B_NMDA + (1. - exp(dt*(( - 1.0 ) / tau_d_NMDA)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau_d_NMDA ) - B_NMDA) ;
   }
  return 0;
}
 
#if NET_RECEIVE_BUFFERING 
#undef t
#define t _nrb_t
static void _net_receive_kernel(double, Point_process*, int _weight_index, double _flag);
static void _net_buf_receive(_NrnThread* _nt) {
  if (!_nt->_ml_list) { return; }
  _Memb_list* _ml = _nt->_ml_list[_mechtype];
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
  /*printf("_net_buf_receive__ProbAMPANMDA_EMS  %d\n", _nt->_id);*/
 
}
 
void _net_receive2 (_NrnThread * _nt, _Memb_list* _ml, int _iml, int _weight_index, double _lflag, double _nrb_t);
void _net_receive (Point_process* _pnt, int _weight_index, double _lflag) {
  _NrnThread* _nt = nrn_threads + _pnt->_tid;
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
 
void _net_receive2 (_NrnThread * _nt, _Memb_list* _ml, int _iml, int _weight_index, double _lflag);
void _net_receive (Point_process* _pnt, int _weight_index, double _lflag)
#endif

{   _Memb_list* _ml;  int _iml;

   _NrnThread* _nt;
   int _tid = _pnt->_tid;
   _nt = nrn_threads + _tid;

   _ml = _nt->_ml_list[_pnt->_type];
   _iml = _pnt->_i_instance;
#if NET_RECEIVE_BUFFERING
   _net_receive2(_nt, _ml, _iml, _weight_index, _lflag, _nrb_t);
#else
   _net_receive2(_nt, _ml, _iml, _weight_index, _lflag);
#endif
}

#if NET_RECEIVE_BUFFERING
void _net_receive2 (_NrnThread * _nt, _Memb_list* _ml, int _iml, int _weight_index, double _lflag, double _nrb_t)
#else
void _net_receive2 (_NrnThread * _nt, _Memb_list* _ml, int _iml, int _weight_index, double _lflag)
#endif
{
   double* _p; Datum* _ppvar;  double v; int _cntml_padded, _cntml_actual; double* _args;
   ThreadDatum* _thread;

   _thread = (ThreadDatum*)0;
   double *_weights = _nt->_weights;
   _args = _weights + _weight_index;
   _cntml_actual = _ml->_nodecount;
   _cntml_padded = _ml->_nodecount_padded;

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
 _tsav = t; {
   double _lresult ;
 _args[1] = _args[0] ;
   _args[2] = _args[0] * NMDA_ratio ;
   if (  ! ( _args[0] > 0.0 ) ) {
     
/*VERBATIM*/
        return;
 }
   if ( Fac > 0.0 ) {
     u = u * exp ( - ( t - tsyn_fac ) / Fac ) ;
     }
   else {
     u = Use ;
     }
   if ( Fac > 0.0 ) {
     u = u + Use * ( 1.0 - u ) ;
     }
   tsyn_fac = t ;
   if ( Rstate  == 0.0 ) {
     _args[3] = exp ( - ( t - tsyn ) / Dep ) ;
     _lresult = urand ( _threadargs_ ) ;
     if ( _lresult > _args[3] ) {
       Rstate = 1.0 ;
       if ( verboseLevel > 0.0 ) {
         printf ( "Recovered! %f at time %g: Psurv = %g, urand=%g\n" , synapseID , t , _args[3] , _lresult ) ;
         }
       }
     else {
       tsyn = t ;
       if ( verboseLevel > 0.0 ) {
         printf ( "Failed to recover! %f at time %g: Psurv = %g, urand=%g\n" , synapseID , t , _args[3] , _lresult ) ;
         }
       }
     }
   if ( Rstate  == 1.0 ) {
     _lresult = urand ( _threadargs_ ) ;
     if ( _lresult < u ) {
       tsyn = t ;
       Rstate = 0.0 ;
       A_AMPA = A_AMPA + _args[1] * factor_AMPA ;
       B_AMPA = B_AMPA + _args[1] * factor_AMPA ;
       A_NMDA = A_NMDA + _args[2] * factor_NMDA ;
       B_NMDA = B_NMDA + _args[2] * factor_NMDA ;
       if ( verboseLevel > 0.0 ) {
         printf ( "Release! %f at time %g: vals %g %g %g %g\n" , synapseID , t , A_AMPA , _args[1] , factor_AMPA , _args[0] ) ;
         }
       }
     else {
       if ( verboseLevel > 0.0 ) {
         printf ( "Failure! %f at time %g: urand = %g\n" , synapseID , t , _lresult ) ;
         }
       }
     }
   } 
#if NET_RECEIVE_BUFFERING
#undef t
#define t _nt->_t
#endif
 }
 
static void _net_init(Point_process* _pnt, int _weight_index, double _lflag) {
   
   double* _p; Datum* _ppvar; ThreadDatum* _thread; 
   _Memb_list* _ml; int _cntml_padded, _cntml_actual; int _iml; double* _args;
 
   _NrnThread* _nt;
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
 }
 
static int  setRNG ( _threadargsproto_ ) {
   
/*VERBATIM*/
    {
#if !NRNBBCORE
	nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
	if (*pv) {
		nrnran123_deletestream(*pv);
		*pv = (nrnran123_State*)0;
	} 
	if (ifarg(2)) {
		*pv = nrnran123_newstream((uint32_t)*getarg(1), (uint32_t)*getarg(2));
	}
#endif
    }
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_setRNG(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 setRNG ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
double urand ( _threadargsproto_ ) {
   double _lurand;
 
/*VERBATIM*/
	double value;
	if (_p_rng) {
		/*
		:Supports separate independent but reproducible streams for
		: each instance.
		*/
		value = nrnran123_negexp((nrnran123_State*)_p_rng);
		//printf("random stream for this simulation = %lf\n",value);
		return value;
	}else{
	        value = 0.0;
//		assert(0);
	}
	_lurand = value;
 
return _lurand;
 }
 
#if 0 /*BBCORE*/
 
static double _hoc_urand(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  urand ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
/*VERBATIM*/
static void bbcore_write(double* x, int* d, int* xx, int* offset, _threadargsproto_) {
	if (d) {
		uint32_t* di = ((uint32_t*)d) + *offset;
		nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
		nrnran123_getids(*pv, di, di+1);
//printf("ProbAMPANMDA_EMS bbcore_write %d %d\n", di[0], di[1]);
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
//printf("ProbAMPANMDA_EMS bbcore_read %d %d\n", di[0], di[1]);
	*offset += 2;
}
 
double toggleVerbose ( _threadargsproto_ ) {
   double _ltoggleVerbose;
 verboseLevel = 1.0 - verboseLevel ;
   
return _ltoggleVerbose;
 }
 
#if 0 /*BBCORE*/
 
static double _hoc_toggleVerbose(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  toggleVerbose ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/

static void initmodel(_threadargsproto_) {
  int _i; double _save;{
  A_NMDA = A_NMDA0;
  A_AMPA = A_AMPA0;
  B_NMDA = B_NMDA0;
  B_AMPA = B_AMPA0;
 {
   double _ltp_AMPA , _ltp_NMDA ;
 
/*VERBATIM*/
          if (_p_rng)
          {
            nrnran123_setseq((nrnran123_State*)_p_rng, 0, 0);
          }
 Rstate = 1.0 ;
   tsyn_fac = 0.0 ;
   tsyn = 0.0 ;
   u = u0 ;
   A_AMPA = 0.0 ;
   B_AMPA = 0.0 ;
   A_NMDA = 0.0 ;
   B_NMDA = 0.0 ;
   _ltp_AMPA = ( tau_r_AMPA * tau_d_AMPA ) / ( tau_d_AMPA - tau_r_AMPA ) * log ( tau_d_AMPA / tau_r_AMPA ) ;
   _ltp_NMDA = ( tau_r_NMDA * tau_d_NMDA ) / ( tau_d_NMDA - tau_r_NMDA ) * log ( tau_d_NMDA / tau_r_NMDA ) ;
   factor_AMPA = - exp ( - _ltp_AMPA / tau_r_AMPA ) + exp ( - _ltp_AMPA / tau_d_AMPA ) ;
   factor_AMPA = 1.0 / factor_AMPA ;
   factor_NMDA = - exp ( - _ltp_NMDA / tau_r_NMDA ) + exp ( - _ltp_NMDA / tau_d_NMDA ) ;
   factor_NMDA = 1.0 / factor_NMDA ;
   A_AMPA_step = exp ( dt * ( ( - 1.0 ) / tau_r_AMPA ) ) ;
   B_AMPA_step = exp ( dt * ( ( - 1.0 ) / tau_d_AMPA ) ) ;
   A_NMDA_step = exp ( dt * ( ( - 1.0 ) / tau_r_NMDA ) ) ;
   B_NMDA_step = exp ( dt * ( ( - 1.0 ) / tau_d_NMDA ) ) ;
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
 _tsav = -1e20;
    _v = _vec_v[_nd_idx];
    _PRCELLSTATE_V
 v = _v;
 _PRCELLSTATE_V
 initmodel(_threadargs_);
}
}

static double _nrn_current(_threadargsproto_, double _v){double _current=0.;v=_v;{ {
   mggate = 1.0 / ( 1.0 + exp ( 0.062 * - ( v ) ) * ( mg / 3.57 ) ) ;
   g_AMPA = gmax * ( B_AMPA - A_AMPA ) ;
   g_NMDA = gmax * ( B_NMDA - A_NMDA ) * mggate ;
   g = g_AMPA + g_NMDA ;
   i_AMPA = g_AMPA * ( v - e ) ;
   i_NMDA = g_NMDA * ( v - e ) ;
   i = i_AMPA + i_NMDA ;
   }
 _current += i;

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
double * _vec_d = _nt->_actual_d;
double * _vec_shadow_rhs = _nt->_shadow_rhs;
double * _vec_shadow_d = _nt->_shadow_d;

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
 {   state(_threadargs_);
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
 
 _slist1 = (int*)malloc(sizeof(int)*4);
 _dlist1 = (int*)malloc(sizeof(int)*4);
 _slist1[0] = &(A_AMPA) - _p;  _dlist1[0] = &(DA_AMPA) - _p;
 _slist1[1] = &(B_AMPA) - _p;  _dlist1[1] = &(DB_AMPA) - _p;
 _slist1[2] = &(A_NMDA) - _p;  _dlist1[2] = &(DA_NMDA) - _p;
 _slist1[3] = &(B_NMDA) - _p;  _dlist1[3] = &(DB_NMDA) - _p;
 #pragma acc enter data copyin(_slist1[0:4])
 #pragma acc enter data copyin(_dlist1[0:4])

_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif
