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
 
#define nrn_init _nrn_init__BinReportHelper
#define nrn_cur _nrn_cur__BinReportHelper
#define _nrn_current _nrn_current__BinReportHelper
#define nrn_jacob _nrn_jacob__BinReportHelper
#define nrn_state _nrn_state__BinReportHelper
#define initmodel initmodel__BinReportHelper
#define _net_receive _net_receive__BinReportHelper
#define nrn_state_launcher nrn_state_BinReportHelper_launcher
#define nrn_cur_launcher nrn_cur_BinReportHelper_launcher
#define nrn_jacob_launcher nrn_jacob_BinReportHelper_launcher 
#define clear clear_BinReportHelper 
#define disable_auto_flush disable_auto_flush_BinReportHelper 
#define flush flush_BinReportHelper 
#define make_comm make_comm_BinReportHelper 
#define pre_savestate pre_savestate_BinReportHelper 
#define restorestate restorestate_BinReportHelper 
#define restoretime restoretime_BinReportHelper 
#define savestate savestate_BinReportHelper 
#define set_steps_to_buffer set_steps_to_buffer_BinReportHelper 
 
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
#define activeStep _p[0*_STRIDE]
#define initialStep _p[1*_STRIDE]
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
 static double _hoc_clear();
 static double _hoc_disable_auto_flush();
 static double _hoc_flush();
 static double _hoc_make_comm();
 static double _hoc_pre_savestate();
 static double _hoc_redirect();
 static double _hoc_restorestate();
 static double _hoc_restoretime();
 static double _hoc_savestate();
 static double _hoc_set_steps_to_buffer();
 
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
 "clear", _hoc_clear,
 "disable_auto_flush", _hoc_disable_auto_flush,
 "flush", _hoc_flush,
 "make_comm", _hoc_make_comm,
 "pre_savestate", _hoc_pre_savestate,
 "redirect", _hoc_redirect,
 "restorestate", _hoc_restorestate,
 "restoretime", _hoc_restoretime,
 "savestate", _hoc_savestate,
 "set_steps_to_buffer", _hoc_set_steps_to_buffer,
 0, 0
};
 
#endif /*BBCORE*/
#define redirect redirect_BinReportHelper
 inline double redirect( _threadargsproto_ );
 /* declare global and static user variables */
#define Dt Dt_BinReportHelper
 double Dt = 0.1;
 #pragma acc declare copyin (Dt)
 
static void _acc_globals_update() {
 #pragma acc update device (Dt) if(nrn_threads->compute_gpu)
 }
 
#if 0 /*BBCORE*/
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "Dt_BinReportHelper", "ms",
 0,0
};
 
#endif /*BBCORE*/
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "Dt_BinReportHelper", &Dt_BinReportHelper,
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
"BinReportHelper",
 "activeStep",
 "initialStep",
 0,
 0,
 0,
 0};
 
static void nrn_alloc(double* _p, Datum* _ppvar, int _type) {
 
#if 0 /*BBCORE*/
 	/*initialize range parameters*/
 	activeStep = 0;
 	initialStep = 0;
 
#endif /* BBCORE */
 
#if 0 /*BBCORE*/
if (!nrn_point_prop_) {_constructor(_prop);}
#endif
 
}
 static void _initlists();
 
#define _tqitem &(_nt->_vdata[_ppvar[2*_STRIDE]])
 static void _net_receive(Point_process*, int, double);
 
#define _psize 4
#define _ppsize 3
 void _BinReportHelper_reg() {
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
  hoc_register_dparam_semantics(_mechtype, 2, "netsend");
 add_nrn_artcell(_mechtype, 2);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, NULL);
 }
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static inline int clear(_threadargsproto_);
static inline int disable_auto_flush(_threadargsproto_);
static inline int flush(_threadargsproto_);
static inline int make_comm(_threadargsproto_);
static inline int pre_savestate(_threadargsproto_);
static inline int restorestate(_threadargsproto_);
static inline int restoretime(_threadargsproto_);
static inline int savestate(_threadargsproto_);
static inline int set_steps_to_buffer(_threadargsproto_);
 } using namespace coreneuron; 
/*VERBATIM*/
#ifndef CORENEURON_BUILD
#ifndef DISABLE_REPORTINGLIB
#include "reportinglib/Records.h"
#include "reportinglib/isc/iscAPI.h"
#include <mpi.h>

extern double* hoc_pgetarg(int iarg);
extern double* getarg(int iarg);
extern char* gargstr(int iarg);
extern int hoc_is_str_arg(int iarg);
extern int ifarg(int iarg);
extern double chkarg(int iarg, double low, double high);
extern double* nrn_recalc_ptr(double*);
extern void nrn_register_recalc_ptr_callback(void (*f)(void));

extern Object** hoc_objgetarg();
extern void* bbss_buffer_counts( int*, int**, int**, int* );
extern void bbss_save_global( void*, char*, int );
extern void bbss_restore_global( void*, char*, int );
extern void bbss_save( void*, int, char*, int );
extern void bbss_restore( void*, int, int, char*, int );
extern void bbss_save_done( void* );
extern void bbss_restore_done( void* );

extern int nrnmpi_myid;

void refreshPointers() { //callback function to update data locations before runtime
	records_refresh_pointers(nrn_recalc_ptr); //tell bin report library to update its pointers using nrn_recalc_ptr function
        isc_refresh_pointers(nrn_recalc_ptr);
}
#endif

//int len=0, *gids=0, *sizes=0, global_size=0, pieceCount=0;
void *bbss_ref = NULL;

#endif
 namespace coreneuron { 
static void _net_receive (Point_process* _pnt, int _weight_index, double _lflag) 
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
#ifndef CORENEURON_BUILD
#ifndef DISABLE_REPORTINGLIB
	records_rec(activeStep);
        isc_send_data(activeStep);

	activeStep++;
#endif
#endif
 artcell_net_send ( _tqitem, _weight_index, _pnt, t +  Dt , 1.0 ) ;
   } 
#if NET_RECEIVE_BUFFERING
#undef t
#define t _nt->_t
#endif
 }
 
static int  make_comm ( _threadargsproto_ ) {
   
/*VERBATIM*/
{
#ifndef CORENEURON_BUILD
#ifndef DISABLE_REPORTINGLIB
	records_setup_communicator();
#endif
#endif
}
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_make_comm(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 make_comm ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
static int  disable_auto_flush ( _threadargsproto_ ) {
   
/*VERBATIM*/
#ifndef CORENEURON_BUILD
#ifndef DISABLE_REPORTINGLIB
    records_set_auto_flush(0);
#endif
#endif
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_disable_auto_flush(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 disable_auto_flush ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
static int  set_steps_to_buffer ( _threadargsproto_ ) {
   
/*VERBATIM*/
#ifndef CORENEURON_BUILD
#ifndef DISABLE_REPORTINGLIB
    int nsteps = (int) *getarg(1);
    records_set_steps_to_buffer( nsteps );
#endif
#endif
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_set_steps_to_buffer(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 set_steps_to_buffer ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
static int  flush ( _threadargsproto_ ) {
   
/*VERBATIM*/
#ifndef CORENEURON_BUILD
#ifndef DISABLE_REPORTINGLIB
        // Note: flush uses actual time (t) whereas recData uses timestep.  Should try to only use one or the other in the future
	records_flush( t );
#endif
#endif
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_flush(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 flush ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
static int  pre_savestate ( _threadargsproto_ ) {
   
/*VERBATIM*/
{
#ifndef CORENEURON_BUILD
#ifndef DISABLE_REPORTINGLIB
    int len, *gids, *sizes, global_size;
    char* gbuffer;
    char* saveFile = gargstr(1);

    //get sizes
    bbss_ref = bbss_buffer_counts( &len, &gids, &sizes, &global_size );

    //pass arrays to bin report library for header creation
    gbuffer = records_saveinit( saveFile, len, gids, sizes, global_size );

    //have neuron fill in global data for rank 0
    if( global_size ) {
        bbss_save_global( bbss_ref, gbuffer, global_size );
        records_saveglobal();
    }

    //for each gid, get the buffer from the bin report lib and give to NEURON layer
    int gidIndex=0;
    for( gidIndex=0; gidIndex<len; gidIndex++ ) {
        char *buffer = records_savebuffer( gids[gidIndex] );
        bbss_save( bbss_ref, gids[gidIndex], buffer, sizes[gidIndex] );
    }

    if(len) {
        free(gids);
        free(sizes);
    }

#endif
#endif
}
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_pre_savestate(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 pre_savestate ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
static int  savestate ( _threadargsproto_ ) {
   
/*VERBATIM*/
#ifndef CORENEURON_BUILD
#ifndef DISABLE_REPORTINGLIB

    if(nrnmpi_myid == 0) {
        printf(" Call to ReportingLib for MPI I/O\n");
    }

    //all buffers populated, have lib execute final MPI-IO operations
    records_savestate();

    //clean up -> I need to free some space.  If they were alloced using 'new', then need report lib to do it
    bbss_save_done( bbss_ref );
#endif
#endif
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_savestate(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 savestate ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
static int  restoretime ( _threadargsproto_ ) {
   
/*VERBATIM*/
#ifndef CORENEURON_BUILD
#ifndef DISABLE_REPORTINGLIB
    void *bbss_ref = NULL;
    int len=0, *gids, *sizes, global_size, pieceCount;
    char* gbuffer = NULL;
    char* saveFile = gargstr(1);

    //get sizes - actually care about gid info for restore
    if( len == 0 ) {
        bbss_ref = bbss_buffer_counts( &len, &gids, &sizes, &global_size );
    }

    // initialize counts, offsets, and get global data - all cpus must load global data unlike with save
    gbuffer = records_restoreinit( saveFile, &global_size );
    bbss_restore_global( bbss_ref, gbuffer, global_size );

    if(len) {
        free(gids);
        free(sizes);
    }
#endif
#endif
 initialStep = t / Dt ;
    return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_restoretime(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 restoretime ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
static int  restorestate ( _threadargsproto_ ) {
   
/*VERBATIM*/
#ifndef CORENEURON_BUILD
#ifndef DISABLE_REPORTINGLIB
{
    void *bbss_ref = NULL;
    int len=0, *gids, *sizes, global_size, pieceCount;
    char* gbuffer = NULL;
    char *saveFile = gargstr(1);

    //get sizes - actually care about gid info for restore
    if( len == 0 ) {
        bbss_ref = bbss_buffer_counts( &len, &gids, &sizes, &global_size );
    }

    // initialize counts, offsets, and get global data - all cpus must load global data unlike with save
    gbuffer = records_restoreinit( saveFile, &global_size );
    bbss_restore_global( bbss_ref, gbuffer, global_size );

    int nbytes = 0, gidIndex=0;
    //for each gid, get the buffer from the bin report lib and give to NEURON layer
    for( gidIndex=0; gidIndex<len; gidIndex++ ) {
        if( gids[gidIndex] != 0 ) {
            gbuffer = records_restore( gids[gidIndex], &pieceCount, &nbytes );
            //printf( "restore %d with %d pieces in %d bytes\n", gids[gidIndex], pieceCount, nbytes );
            bbss_restore( bbss_ref, gids[gidIndex], pieceCount, gbuffer, nbytes );
        }
    }

    //clean up -> I need to free some space.  If they were alloced using 'new', then need report lib to do it
    bbss_restore_done( bbss_ref );

    if(len) {
        free(gids);
        free(sizes);
    }
}
#endif
#endif
 activeStep = t / Dt ;
    return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_restorestate(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 restorestate ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
double redirect ( _threadargsproto_ ) {
   double _lredirect;
 
/*VERBATIM*/
{
#ifndef CORENEURON_BUILD
#ifndef DISABLE_REPORTINGLIB
    FILE *fout;
    char fname[128];

    int mpi_size, mpi_rank;

    // get MPI info
    MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);

    sprintf( fname, "NodeFiles/%d.%dnode.out", mpi_rank, mpi_size );
    fout = freopen( fname, "w", stdout );
    if( !fout ) {
        fprintf( stderr, "failed to redirect.  Terminating\n" );
        exit(0);
    }

    sprintf( fname, "NodeFiles/%d.%dnode.err", mpi_rank, mpi_size );
    fout = freopen( fname, "w", stderr );
    setbuf( fout, NULL );
#endif
#endif
}
 
return _lredirect;
 }
 
#if 0 /*BBCORE*/
 
static double _hoc_redirect(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  redirect ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
static int  clear ( _threadargsproto_ ) {
   
/*VERBATIM*/
#ifndef CORENEURON_BUILD
#ifndef DISABLE_REPORTINGLIB
        records_clear();
#endif
#endif
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_clear(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 clear ( _threadargs_ );
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
{
#ifndef CORENEURON_BUILD
#ifndef DISABLE_REPORTINGLIB
    if(ifarg(1)) {
        Dt = *getarg(1);

        records_set_atomic_step(Dt);
        isc_set_sim_dt(Dt);

        nrn_register_recalc_ptr_callback( refreshPointers );
    }
#endif
#endif
}
 }
 
}
}
#endif /*BBCORE*/

static inline void initmodel(_threadargsproto_) {
  int _i; double _save;{
 {
   activeStep = initialStep ;
   artcell_net_send ( _tqitem, -1, (Point_process*) _nt->_vdata[_ppvar[1*_STRIDE]], t +  initialStep * Dt , 1.0 ) ;
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
 _tsav = -1e20;
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
