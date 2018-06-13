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
 
#define _thread_present_ /**/ , _slist1[0:1], _dlist1[0:1] 
 
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
 
#define nrn_init _nrn_init__ExpSyn
#define nrn_cur _nrn_cur__ExpSyn
#define nrn_cur_parallel _nrn_cur_parallel__ExpSyn
#define _nrn_current _nrn_current__ExpSyn
#define nrn_jacob _nrn_jacob__ExpSyn
#define nrn_state _nrn_state__ExpSyn
#define initmodel initmodel__ExpSyn
#define _net_receive _net_receive__ExpSyn
#define _net_receive2 _net_receive2__ExpSyn
#define nrn_state_launcher nrn_state_ExpSyn_launcher
#define nrn_cur_launcher nrn_cur_ExpSyn_launcher
#define nrn_jacob_launcher nrn_jacob_ExpSyn_launcher 
#if NET_RECEIVE_BUFFERING
#define _net_buf_receive _net_buf_receive_ExpSyn
static void _net_buf_receive(_NrnThread*);
#endif
 
#define state state_ExpSyn
#define _ode_matsol1 _nrn_ode_matsol1__ExpSyn
#define _ode_spec1 _nrn_ode_spec1__ExpSyn
 
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
#define tau _p[0*_STRIDE]
#define e _p[1*_STRIDE]
#define i _p[2*_STRIDE]
#define g _p[3*_STRIDE]
#define Dg _p[4*_STRIDE]
#define _v_unused _p[5*_STRIDE]
#define _g_unused _p[6*_STRIDE]
#define _tsav _p[7*_STRIDE]
 
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
 0, 0
};
 
#endif /*BBCORE*/
 /* declare global and static user variables */
 
static void _acc_globals_update() {
 }
 
#if 0 /*BBCORE*/
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "tau", 1e-09, 1e+09,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tau", "ms",
 "e", "mV",
 "g", "uS",
 "i", "nA",
 0,0
};
 
#endif /*BBCORE*/
 static double delta_t = 0.01;
 static double g0 = 0;
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
 
#if 0 /*BBCORE*/
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 
#endif /*BBCORE*/
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "6.2.0",
"ExpSyn",
 "tau",
 "e",
 0,
 "i",
 0,
 "g",
 0,
 0};
 
 void _nrn_ode_state_vars__ExpSyn(short * count, short** var_offsets, short ** dv_offsets)
 {
     *count = 1;
     (*var_offsets) = (short*) malloc(sizeof(short)* *count);
     (*dv_offsets) = (short*) malloc(sizeof(short)* *count);
     (*var_offsets)[0] = 3;
     (*dv_offsets)[0] = 4;
 }

static void nrn_alloc(double* _p, Datum* _ppvar, int _type) {
 
#if 0 /*BBCORE*/
 	/*initialize range parameters*/
 	tau = 0.1;
 	e = 0;
 
#endif /* BBCORE */
 
}
 static void _initlists();
 static void _net_receive(Point_process*, int, double);
 
#define _psize 8
#define _ppsize 2
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*f)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(_threadargsproto_, int));
extern void _cvode_abstol( Symbol**, double*, int);

 void _expsyn_reg() {
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
  hoc_register_prop_size(_mechtype, _psize, _ppsize);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
 
#if NET_RECEIVE_BUFFERING
  hoc_register_net_receive_buffering(_net_buf_receive, _mechtype);
#endif
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, NULL);
 }
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
int _ode_spec1(_threadargsproto_);
/*int _ode_matsol1(_threadargsproto_);*/
 
#define _slist1 _slist1_ExpSyn
int* _slist1;
#pragma acc declare create(_slist1)

#define _dlist1 _dlist1_ExpSyn
int* _dlist1;
#pragma acc declare create(_dlist1)
 static inline int state(_threadargsproto_);
 
/*CVODE*/
int _ode_spec1 (_threadargsproto_) {int _reset = 0; {
   Dg = - g / tau ;
   }
 return _reset;
}

int _ode_matsol1 (_threadargsproto_) {
 Dg = Dg  / (1. - dt*( ( - 1.0 ) / tau )) ;
 return 0;
}
 /*END CVODE*/
 static int state (_threadargsproto_) { {
    g = g + (1. - exp(dt*(( - 1.0 ) / tau)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau ) - g) ;
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
  /*printf("_net_buf_receive__ExpSyn  %d\n", _nt->_id);*/
 
}
 
void _net_receive2 (_NrnThread * _nt, _Memb_list* _ml, int _iml, int _weight_index, double _lflag, double _nrb_t);
static void _net_receive (Point_process* _pnt, int _weight_index, double _lflag) {
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
   double *_weights = _nt->_weights;
   ThreadDatum* _thread;

   _thread = (ThreadDatum*)0;
   _args = _weights + _weight_index;
   _cntml_actual = _ml->_nodecount;
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
   g = g + _args[0] ;
   } 
#if NET_RECEIVE_BUFFERING
#undef t
#define t _nt->_t
#endif
 }

static void initmodel(_threadargsproto_) {
  int _i; double _save;{
  g = g0;
 {
   g = 0.0 ;
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
   i = g * ( v - e ) ;
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
double * _vec_shadow_rhs = _ml->_shadow_rhs;
double * _vec_shadow_d = _ml->_shadow_d;

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
 
 _slist1 = (int*)malloc(sizeof(int)*1);
 _dlist1 = (int*)malloc(sizeof(int)*1);
 _slist1[0] = &(g) - _p;  _dlist1[0] = &(Dg) - _p;
 #pragma acc enter data copyin(_slist1[0:1])
 #pragma acc enter data copyin(_dlist1[0:1])

_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif
