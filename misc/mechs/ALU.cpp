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
 
#define _thread_present_ /**/ 
 
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
 
#define nrn_init _nrn_init__ALU
#define nrn_cur _nrn_cur__ALU
#define _nrn_current _nrn_current__ALU
#define nrn_jacob _nrn_jacob__ALU
#define nrn_state _nrn_state__ALU
#define initmodel initmodel__ALU
#define _net_receive _net_receive__ALU
#define nrn_state_launcher nrn_state_ALU_launcher
#define nrn_cur_launcher nrn_cur_ALU_launcher
#define nrn_jacob_launcher nrn_jacob_ALU_launcher 
#if NET_RECEIVE_BUFFERING
#define _net_buf_receive _net_buf_receive_ALU
static void _net_buf_receive(NrnThread*);
#endif
 
#define average average_ALU 
#define addvar addvar_ALU 
#define constant constant_ALU 
#define setop setop_ALU 
#define summation summation_ALU 
 
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
#define Dt _p[0*_STRIDE]
#define output _p[1*_STRIDE]
#define _v_unused _p[2*_STRIDE]
#define _tsav _p[3*_STRIDE]
 
#ifndef NRN_PRCELLSTATE
#define NRN_PRCELLSTATE 0
#endif
#if NRN_PRCELLSTATE
#define _PRCELLSTATE_V _v_unused = _v;
#define _PRCELLSTATE_G /**/
#else
#define _PRCELLSTATE_V /**/
#define _PRCELLSTATE_G /**/
#endif
#define _nd_area  _nt_data[_ppvar[0*_STRIDE]]
#define _p_ptr	_nt->_vdata[_ppvar[2*_STRIDE]]
 
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
 /* external NEURON variables */
 
#if 0 /*BBCORE*/
 /* declaration of user functions */
 static double _hoc_average();
 static double _hoc_addvar();
 static double _hoc_constant();
 static double _hoc_setop();
 static double _hoc_summation();
 
#endif /*BBCORE*/
 
#define _mechtype _mechtype_ALU
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
 "average", _hoc_average,
 "addvar", _hoc_addvar,
 "constant", _hoc_constant,
 "setop", _hoc_setop,
 "summation", _hoc_summation,
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
 "Dt", "ms",
 0,0
};
 
#endif /*BBCORE*/
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
 
#if 0 /*BBCORE*/
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 
#endif /*BBCORE*/
 
#if 0 /*BBCORE*/
static void _destructor(Prop*);
#endif
 
#if 0 /*BBCORE*/
static void _constructor(Prop*);
#endif
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "6.2.0",
"ALU",
 "Dt",
 "output",
 0,
 0,
 0,
 "ptr",
 0};
 
static void nrn_alloc(double* _p, Datum* _ppvar, int _type) {
 
#if 0 /*BBCORE*/
 	/*initialize range parameters*/
 	Dt = 0.1;
 	output = 0;
 
#endif /* BBCORE */
 
#if 0 /*BBCORE*/
if (!nrn_point_prop_) {_constructor(_prop);}
#endif
 
}
 static void _initlists();
 
#define _tqitem &(_nt->_vdata[_ppvar[3*_STRIDE]])
 
#if NET_RECEIVE_BUFFERING
#undef _tqitem
#define _tqitem _ppvar[3*_STRIDE]
#endif

 static void _net_receive(Point_process*, int, double);
 
#define _psize 4
#define _ppsize 4
 void _ALU_reg() {
	int _vectorized = 1;
  _initlists();
 _mechtype = nrn_get_mechtype(_mechanism[1]);
 if (_mechtype == -1) return;
 _nrn_layout_reg(_mechtype, LAYOUT);
 
#if 0 /*BBCORE*/
 
#endif /*BBCORE*/
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,NULL, NULL, NULL, nrn_init,
	 hoc_nrnpointerindex,
	 NULL/*_hoc_create_pnt*/, NULL/*_hoc_destroy_pnt*/, /*_member_func,*/
	 1);
     #if 0 /*BBCORE*/
    register_destructor(_destructor);
    #endif
   hoc_reg_bbcore_read(_mechtype, bbcore_read);
   hoc_reg_bbcore_write(_mechtype, bbcore_write);
  hoc_register_prop_size(_mechtype, _psize, _ppsize);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "bbcorepointer");
  hoc_register_dparam_semantics(_mechtype, 3, "netsend");
 
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
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static inline int average(_threadargsproto_);
static inline int addvar(_threadargsproto_);
static inline int constant(_threadargsproto_);
static inline int setop(_threadargsproto_);
static inline int summation(_threadargsproto_);
 } using namespace coreneuron; 
/*VERBATIM*/

#ifndef CORENEURON_BUILD
extern double* hoc_pgetarg(int iarg);
extern double* getarg(int iarg);
extern int ifarg(int iarg);

typedef struct {
    //! list of pointers to hoc variables
	double** ptrs_;

    /*! list of scalars to apply to corresponding variables; useful for making units of variables
     * from different sources consistent (e.g. i current sources may be distributed, mA/cm^2, or point processes, nA)
     */
    double * scalars_;

    //! number of elements stored in the vectors
	int np_;

    //! number of slots allocated to the vectors
	int psize_;

    //! function pointer to execute when net_receive is triggered
    int (*process)(_threadargsproto_);
} Info;

#define INFOCAST Info** ip = (Info**)(&(_p_ptr))

#define dp double*
extern void nrn_register_recalc_ptr_callback(void (*f)());
extern Point_process* ob2pntproc(Object*);
extern double* nrn_recalc_ptr(double*);

static void recalcptr(Info* info, int cnt, double** old_vp, double* new_v) {
        int i;
        /*printf("recalcptr np_=%d %s\n", info->np_, info->path_);*/

}
static void recalc_ptr_callback() {
        Symbol* sym;
        int i;
        hoc_List* instances;
        hoc_Item* q;
        /*printf("ASCIIrecord.mod recalc_ptr_callback\n");*/
        /* hoc has a list of the ASCIIRecord instances */
        sym = hoc_lookup("ALU");
        instances = sym->u.template->olist;
        ITERATE(q, instances) {
                Info* InfoPtr;
                Point_process* pnt;
                Object* o = OBJ(q);
                /*printf("callback for %s\n", hoc_object_name(o));*/
                pnt = ob2pntproc(o);
                Datum* _ppvar = pnt->_prop->dparam;
                INFOCAST;
                InfoPtr = *ip;
                        for (i=0; i < InfoPtr->np_; ++i)
                                InfoPtr->ptrs_[i] =  nrn_recalc_ptr(InfoPtr->ptrs_[i]);

        }
}
#endif
 namespace coreneuron { 
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
  /*printf("_net_buf_receive__ALU  %d\n", _nt->_id);*/
 
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
   
/*VERBATIM*/
{
#ifndef CORENEURON_BUILD
    INFOCAST;
    Info* info = *ip;
	info->process(_threadargs_);
#endif
}
 
#if NET_RECEIVE_BUFFERING
    _net_send_buffering(_ml->_net_send_buffer, 0, _tqitem, _weight_index, _ppvar[1*_STRIDE], t +  Dt , 1.0 );
#else
 net_send ( _tqitem, _weight_index, _pnt, t +  Dt , 1.0 ) ;
   
#endif
 } 
#if NET_RECEIVE_BUFFERING
#undef t
#define t _nt->_t
#endif
 }
 
static int  addvar ( _threadargsproto_ ) {
   
/*VERBATIM*/
{
#ifndef CORENEURON_BUILD
    INFOCAST;
    Info* info = *ip;
	if (info->np_ >= info->psize_) {
		info->psize_ += 10;
		info->ptrs_ = (double**) hoc_Erealloc(info->ptrs_, info->psize_*sizeof(double*)); hoc_malchk();
        info->scalars_ = (double*) hoc_Erealloc(info->scalars_, info->psize_*sizeof(double)); hoc_malchk();
	}

	info->ptrs_[info->np_] = hoc_pgetarg(1);
    if( ifarg(2)) {
        info->scalars_[info->np_] = *getarg(2);
    } else {
        info->scalars_[info->np_] = 1;
    }

	++info->np_;
    //printf("I have %d values.. (new = %g * %g)\n", info->np_, *(info->ptrs_[info->np_-1]), info->scalars_[info->np_-1] );
#endif
}
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_addvar(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 addvar ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
static int  constant ( _threadargsproto_ ) {
   
/*VERBATIM*/
{
#ifndef CORENEURON_BUILD
    INFOCAST;
    Info* info = *ip;
    if( info->np_ > 0 ) {
        output = info->scalars_[0];
    } else {
        output = 0;
    }
#endif
}
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_constant(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 constant ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
static int  average ( _threadargsproto_ ) {
   
/*VERBATIM*/
{
#ifndef CORENEURON_BUILD
    INFOCAST;
    Info* info = *ip;
	int i;
	double n = 0;
	for (i=0; i < info->np_; ++i) {
      //  printf("%f", (*info->ptrs_[i] * info->scalars_[i]) );
		n += (*info->ptrs_[i] * info->scalars_[i]);
	}
    //printf("\n");
//	output = n/info->np_;
	if (info->np_ > 0)
	  output = n/info->np_;
	else output = 0;
#endif
}
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_average(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 average ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
static int  summation ( _threadargsproto_ ) {
   
/*VERBATIM*/
{
#ifndef CORENEURON_BUILD
    INFOCAST; Info* info = *ip;
	int i;
	double n = 0;
	for (i=0; i < info->np_; ++i) {
        //printf("%f = %f * %f\n", (*info->ptrs_[i] * info->scalars_[i]), *info->ptrs_[i], info->scalars_[i] );
		n += (*info->ptrs_[i] * info->scalars_[i]);
	}

    output = n;
#endif
}
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_summation(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 summation ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
static int  setop ( _threadargsproto_ ) {
   
/*VERBATIM*/
{
#ifndef CORENEURON_BUILD
    INFOCAST; Info* info = *ip;

    char *opname = NULL;
    if (!hoc_is_str_arg(1)) {
        exit(0);
    }

    opname = gargstr(1);
    if( strcmp( opname, "summation" ) == 0 ) {
        info->process = &summation;
    } else if ( strcmp( opname, "average" ) == 0 ) {
        info->process = &average;
    } else if ( strcmp( opname, "constant" ) == 0 ) {
        info->process = &constant;
    } else {
        fprintf( stderr, "Error: unknown operation '%s' for ALU object.  Terminating.\n", opname );
        exit(0);
    }
#endif
}
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_setop(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 setop ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 } using namespace coreneuron; 
/*VERBATIM*/
/** not executed in coreneuron and hence need empty stubs only */
static void bbcore_write(double* x, int* d, int* xx, int* offset, _threadargsproto_) {
}
static void bbcore_read(double* x, int* d, int* xx, int* offset, _threadargsproto_) {
}
 namespace coreneuron { 
#if 0 /*BBCORE*/
static _constructor(_prop)
	Prop *_prop; double* _p; Datum* _ppvar; ThreadDatum* _thread;
{
	_thread = (Datum*)0;
	_p = _prop->param; _ppvar = _prop->dparam;
{
 {
   
/*VERBATIM*/
{
#ifndef CORENEURON_BUILD
        static int first = 1;
        if (first) {
                first = 0;
                nrn_register_recalc_ptr_callback(recalc_ptr_callback);
        }

	INFOCAST;
	Info* info = (Info*)hoc_Emalloc(sizeof(Info)); hoc_malchk();
	info->psize_ = 10;
	info->ptrs_ = (double**)hoc_Ecalloc(info->psize_, sizeof(double*)); hoc_malchk();
    info->scalars_ = (double*)hoc_Ecalloc(info->psize_, sizeof(double)); hoc_malchk();
	info->np_ = 0;
	*ip = info;

    if (ifarg(2)) {
         Dt = *getarg(2);
    }

    //default operation is average
    info->process = &average;
#endif
}
 }
 
}
}
#endif /*BBCORE*/
 
#if 0 /*BBCORE*/
static _destructor(_prop)
	Prop *_prop; double* _p; Datum* _ppvar; ThreadDatum* _thread;
{
	_thread = (Datum*)0;
	_p = prop->param; _ppvar = _prop->dparam;
{
 {
   
/*VERBATIM*/
{
#ifndef CORENEURON_BUILD
	INFOCAST;
    Info* info = *ip;
	free(info->ptrs_);
	free(info);
#endif
}
 }
 
}
}
#endif /*BBCORE*/

static inline void initmodel(_threadargsproto_) {
  int _i; double _save;{
  Memb_list* _ml = _nt->_ml_list[_mechtype];
 {
   
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

static double _nrn_current(_threadargsproto_, double _v){double _current=0.;v=_v;{
} return _current;
}

#if defined(ENABLE_CUDA_INTERFACE) && defined(_OPENACC)
  void nrn_state_launcher(NrnThread*, Memb_list*, int, int);
  void nrn_jacob_launcher(NrnThread*, Memb_list*, int, int);
  void nrn_cur_launcher(NrnThread*, Memb_list*, int, int);
#endif


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
}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
 int _cntml_actual=1;
 int _cntml_padded=1;
 int _iml=0;
  if (!_first) return;
_first = 0;
}
} // namespace coreneuron_lib
