/* Created by Language version: 6.2.0 */

#undef DISABLE_OPENACC
#define DISABLE_OPENACC

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
 
#define nrn_init _nrn_init__MemUsage
#define nrn_cur _nrn_cur__MemUsage
#define _nrn_current _nrn_current__MemUsage
#define nrn_jacob _nrn_jacob__MemUsage
#define nrn_state _nrn_state__MemUsage
#define initmodel initmodel__MemUsage
#define _net_receive _net_receive__MemUsage
#define nrn_state_launcher nrn_state_MemUsage_launcher
#define nrn_cur_launcher nrn_cur_MemUsage_launcher
#define nrn_jacob_launcher nrn_jacob_MemUsage_launcher 
#define print_node_mem_usage print_node_mem_usage_MemUsage 
#define print_mem_usage print_mem_usage_MemUsage 
 
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
#define _v_unused _p[0*_STRIDE]
 
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
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 static int hoc_nrnpointerindex =  -1;
 static ThreadDatum* _extcall_thread;
 /* external NEURON variables */
 
#if 0 /*BBCORE*/
 /* declaration of user functions */
 static double _hoc_print_node_mem_usage();
 static double _hoc_print_mem_usage();
 
#endif /*BBCORE*/
 static int _mechtype;
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
 "print_node_mem_usage", _hoc_print_node_mem_usage,
 "print_mem_usage", _hoc_print_mem_usage,
 0, 0
};
 
#endif /*BBCORE*/
 /* declare global and static user variables */
#define avgUsageMB avgUsageMB_MemUsage
 double avgUsageMB = 0;
 #pragma acc declare copyin (avgUsageMB)
#define maxUsageMB maxUsageMB_MemUsage
 double maxUsageMB = 0;
 #pragma acc declare copyin (maxUsageMB)
#define minUsageMB minUsageMB_MemUsage
 double minUsageMB = 0;
 #pragma acc declare copyin (minUsageMB)
#define rank rank_MemUsage
 double rank = 0;
 #pragma acc declare copyin (rank)
#define size size_MemUsage
 double size = 0;
 #pragma acc declare copyin (size)
#define stdevUsageMB stdevUsageMB_MemUsage
 double stdevUsageMB = 0;
 #pragma acc declare copyin (stdevUsageMB)
 
static void _acc_globals_update() {
 #pragma acc update device (avgUsageMB) if(nrn_threads->compute_gpu)
 #pragma acc update device (maxUsageMB) if(nrn_threads->compute_gpu)
 #pragma acc update device (minUsageMB) if(nrn_threads->compute_gpu)
 #pragma acc update device (rank) if(nrn_threads->compute_gpu)
 #pragma acc update device (size) if(nrn_threads->compute_gpu)
 #pragma acc update device (stdevUsageMB) if(nrn_threads->compute_gpu)
 }
 
#if 0 /*BBCORE*/
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 0,0
};
 
#endif /*BBCORE*/
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "minUsageMB_MemUsage", &minUsageMB_MemUsage,
 "maxUsageMB_MemUsage", &maxUsageMB_MemUsage,
 "avgUsageMB_MemUsage", &avgUsageMB_MemUsage,
 "stdevUsageMB_MemUsage", &stdevUsageMB_MemUsage,
 "rank_MemUsage", &rank_MemUsage,
 "size_MemUsage", &size_MemUsage,
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
static void _constructor(Prop*);
#endif
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "6.2.0",
"MemUsage",
 0,
 0,
 0,
 0};
 
static void nrn_alloc(double* _p, Datum* _ppvar, int _type) {
 
#if 0 /*BBCORE*/
 	/*initialize range parameters*/
 
#endif /* BBCORE */
 
#if 0 /*BBCORE*/
if (!nrn_point_prop_) {_constructor(_prop);}
#endif
 
}
 static void _initlists();
 
#define _psize 1
#define _ppsize 2
 void _MemUsage_reg() {
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
  hoc_register_prop_size(_mechtype, _psize, _ppsize);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
 add_nrn_artcell(_mechtype, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, NULL);
 }
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static inline int print_node_mem_usage(_threadargsproto_);
static inline int print_mem_usage(_threadargsproto_);
 } using namespace coreneuron; 
/*VERBATIM*/
// HPE mpi fails due to t macro definition, skip mpi.h include
#ifndef CORENEURON_BUILD
#include <mpi.h>
#endif
int nrn_mallinfo(int);
 namespace coreneuron { 
static int  print_mem_usage ( _threadargsproto_ ) {
   
/*VERBATIM*/
#ifndef CORENEURON_BUILD
/**
 * Gather memory usage statistics for all nodes in the network, printing to the console
 * argument 6 to nrn_mallinfo includes memory mapped files (m.hblkhd + m.arena)
 * argument 1 to nrn_mallinfo returns uordblks which is "total size of memory occupied by chunks handed out by malloc"
 */
    double usageMB = (double) nrn_mallinfo(1) / (double) (1024*1024);

    MPI_Reduce( &usageMB, &minUsageMB, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Reduce( &usageMB, &maxUsageMB, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce( &usageMB, &avgUsageMB, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

    avgUsageMB /= size;

    MPI_Bcast( &avgUsageMB, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    double diffSquared = (usageMB-avgUsageMB)*(usageMB-avgUsageMB);
    MPI_Reduce( &diffSquared, &stdevUsageMB, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    stdevUsageMB = sqrt( stdevUsageMB/size);

    if( rank == 0 ) {
        printf( " memusage according to nrn_mallinfo:\tMax = %lfMB\tMin = %lfMB\tMean(Stdev) = %lfMB(%lfMB)\n", maxUsageMB, minUsageMB, avgUsageMB, stdevUsageMB );
    }
#endif
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_print_mem_usage(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 print_mem_usage ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
static int  print_node_mem_usage ( _threadargsproto_ ) {
   
/*VERBATIM*/
#ifndef CORENEURON_BUILD
    double usageMB = (double) nrn_mallinfo(6) / (double) (1024*1024);
    printf( " memusage node %.0lf according to nrn_mallinfo:\t %lfMB\n", rank, usageMB );
#endif
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_print_node_mem_usage(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 print_node_mem_usage ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
#if 0 /*BBCORE*/
static _constructor(_prop)
	Prop *_prop; double* _p; Datum* _ppvar; ThreadDatum* _thread;
{
	_thread = (Datum*)0;
	_p = _prop->param; _ppvar = _prop->dparam;
{
 {
   
/*VERBATIM*/
    int i_rank, i_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &i_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &i_size);
    rank = i_rank;
    size = i_size;
 }
 
}
}
#endif /*BBCORE*/

static inline void initmodel(_threadargsproto_) {
  int _i; double _save;{
 {
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
_PRAGMA_FOR_INIT_ACC_LOOP_
for (_iml = 0; _iml < _cntml_actual; ++_iml) {
#else /* LAYOUT > 1 */ /*AoSoA*/
#error AoSoA not implemented.
for (;;) { /* help clang-format properly indent */
#endif
 initmodel(_threadargs_);
}
  }
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
