/* Created by Language version: 6.2.0 */

#undef DISABLE_OPENACC
#define DISABLE_OPENACC

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
 
#define nrn_init _nrn_init__InhPoissonStim
#define nrn_cur _nrn_cur__InhPoissonStim
#define _nrn_current _nrn_current__InhPoissonStim
#define nrn_jacob _nrn_jacob__InhPoissonStim
#define nrn_state _nrn_state__InhPoissonStim
#define initmodel initmodel__InhPoissonStim
#define _net_receive _net_receive__InhPoissonStim
#define _net_receive2 _net_receive2__InhPoissonStim
#define nrn_state_launcher nrn_state_InhPoissonStim_launcher
#define nrn_cur_launcher nrn_cur_InhPoissonStim_launcher
#define nrn_jacob_launcher nrn_jacob_InhPoissonStim_launcher 
#define generate_next_event generate_next_event_InhPoissonStim 
#define setRate setRate_InhPoissonStim 
#define setTbins setTbins_InhPoissonStim 
#define setRNGs setRNGs_InhPoissonStim 
#define update_time update_time_InhPoissonStim 
 
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
#define duration _p[0*_STRIDE]
#define rmax _p[1*_STRIDE]
#define index _p[2*_STRIDE]
#define curRate _p[3*_STRIDE]
#define start _p[4*_STRIDE]
#define event _p[5*_STRIDE]
#define _v_unused _p[6*_STRIDE]
#define _tsav _p[7*_STRIDE]
 
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
#define _p_uniform_rng	_nt->_vdata[_ppvar[2*_STRIDE]]
#define _p_exp_rng	_nt->_vdata[_ppvar[3*_STRIDE]]
#define _p_vecRate	_nt->_vdata[_ppvar[4*_STRIDE]]
#define _p_vecTbins	_nt->_vdata[_ppvar[5*_STRIDE]]
 
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
 static double _hoc_erand();
 static double _hoc_generate_next_event();
 static double _hoc_setRate();
 static double _hoc_setTbins();
 static double _hoc_setRNGs();
 static double _hoc_urand();
 static double _hoc_update_time();
 
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
 "erand", _hoc_erand,
 "generate_next_event", _hoc_generate_next_event,
 "setRate", _hoc_setRate,
 "setTbins", _hoc_setTbins,
 "setRNGs", _hoc_setRNGs,
 "urand", _hoc_urand,
 "update_time", _hoc_update_time,
 0, 0
};
 
#endif /*BBCORE*/
#define erand erand_InhPoissonStim
#define urand urand_InhPoissonStim
 inline double erand( _threadargsproto_ );
 inline double urand( _threadargsproto_ );
 /* declare global and static user variables */
#define interval_min interval_min_InhPoissonStim
 double interval_min = 1;
 #pragma acc declare copyin (interval_min)
 
static void _acc_globals_update() {
 #pragma acc update device (interval_min) if(nrn_threads->compute_gpu)
 }
 
#if 0 /*BBCORE*/
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "duration", 0, 1e+09,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "duration", "ms",
 0,0
};
 
#endif /*BBCORE*/
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "interval_min_InhPoissonStim", &interval_min_InhPoissonStim,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(double*, Datum*, int);
void nrn_init(_NrnThread*, _Memb_list*, int);
void nrn_state(_NrnThread*, _Memb_list*, int);
 
#if 0 /*BBCORE*/
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 
#endif /*BBCORE*/
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "6.2.0",
"InhPoissonStim",
 "duration",
 0,
 "rmax",
 0,
 0,
 "uniform_rng",
 "exp_rng",
 "vecRate",
 "vecTbins",
 0};
 
static void nrn_alloc(double* _p, Datum* _ppvar, int _type) {
 
#if 0 /*BBCORE*/
 	/*initialize range parameters*/
 	duration = 1e+06;
 
#endif /* BBCORE */
 
}
 static void _initlists();
 
#define _tqitem &(_nt->_vdata[_ppvar[6*_STRIDE]])
 static void _net_receive(Point_process*, int, double);
 
#define _psize 8
#define _ppsize 7
 static void bbcore_read(double *, int*, int*, int*, _threadargsproto_);
 extern void hoc_reg_bbcore_read(int, void(*)(double *, int*, int*, int*, _threadargsproto_));
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*f)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(_threadargsproto_, int));
extern void _cvode_abstol( Symbol**, double*, int);

 void _netstim_inhpoisson_reg() {
	int _vectorized = 1;
  _initlists();
 _mechtype = nrn_get_mechtype(_mechanism[1]);
 if (_mechtype == -1) return;
 _nrn_layout_reg(_mechtype, LAYOUT);
 
#if 0 /*BBCORE*/
 
#endif /*BBCORE*/
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,(void*)0, (void*)0, (void*)0, nrn_init,
	 hoc_nrnpointerindex,
	 NULL/*_hoc_create_pnt*/, NULL/*_hoc_destroy_pnt*/, /*_member_func,*/
	 1);
   hoc_reg_bbcore_read(_mechtype, bbcore_read);
  hoc_register_prop_size(_mechtype, _psize, _ppsize);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "bbcorepointer");
  hoc_register_dparam_semantics(_mechtype, 3, "bbcorepointer");
  hoc_register_dparam_semantics(_mechtype, 4, "bbcorepointer");
  hoc_register_dparam_semantics(_mechtype, 5, "bbcorepointer");
  hoc_register_dparam_semantics(_mechtype, 6, "netsend");
 add_nrn_artcell(_mechtype, 6);
 add_nrn_has_net_event(_mechtype);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, NULL);
 }
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int generate_next_event(_threadargsproto_);
static int setRate(_threadargsproto_);
static int setTbins(_threadargsproto_);
static int setRNGs(_threadargsproto_);
static int update_time(_threadargsproto_);
 
/*VERBATIM*/
extern int ifarg(int iarg);
extern double* vector_vec(void* vv);
extern void* vector_new1(int _i);
extern int vector_capacity(void* vv);
extern void* vector_arg(int iarg);
#if !defined(CORENEURON_BUILD)
double nrn_random_pick(void* r);
#endif
void* nrn_random_arg(int argpos);
 
/*VERBATIM*/
#include "nrnran123.h"
 
static int  generate_next_event ( _threadargsproto_ ) {
   event = 1000.0 / rmax * erand ( _threadargs_ ) ;
   if ( event < 0.0 ) {
     event = 0.0 ;
     }
    return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_generate_next_event(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 generate_next_event ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
static int  setRNGs ( _threadargsproto_ ) {
   
/*VERBATIM*/
{
#if !NRNBBCORE
        nrnran123_State** pv = (nrnran123_State**)(&_p_exp_rng);
        if (*pv) {
                nrnran123_deletestream(*pv);
                *pv = (nrnran123_State*)0;
        } 
        if (ifarg(2)) {
                *pv = nrnran123_newstream((uint32_t)*getarg(1), (uint32_t)*getarg(2));
        }

        pv = (nrnran123_State**)(&_p_uniform_rng);
        if (*pv) {
                nrnran123_deletestream(*pv);
                *pv = (nrnran123_State*)0;
        } 
        if (ifarg(4)) {
                *pv = nrnran123_newstream((uint32_t)*getarg(3), (uint32_t)*getarg(4));
        }
#endif
}
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_setRNGs(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 setRNGs ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
double urand ( _threadargsproto_ ) {
   double _lurand;
 
/*VERBATIM*/
	if (_p_uniform_rng) {
		/*
		:Supports separate independent but reproducible streams for
		: each instance. However, the corresponding hoc Random
		: distribution MUST be set to Random.uniform(0,1)
		*/
		_lurand = nrnran123_dblpick((nrnran123_State*)_p_uniform_rng);
	}else{
          _lurand = 0.0;
  	  hoc_execerror("multithread random in NetStim"," only via hoc Random");
	}
 
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
 
double erand ( _threadargsproto_ ) {
   double _lerand;
 
/*VERBATIM*/
	if (_p_exp_rng) {
		/*
		:Supports separate independent but reproducible streams for
		: each instance. However, the corresponding hoc Random
		: distribution MUST be set to Random.negexp(1)
		*/
		_lerand = nrnran123_negexp((nrnran123_State*)_p_exp_rng);
	}else{
          _lerand = 0.0;
  	  hoc_execerror("multithread random in NetStim"," only via hoc Random");
	}
 
return _lerand;
 }
 
#if 0 /*BBCORE*/
 
static double _hoc_erand(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  erand ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
/*VERBATIM*/
static void bbcore_write(double* dArray, int* iArray, int* doffset, int* ioffset, _threadargsproto_) {
        uint32_t dsize = 0;
        if (_p_vecRate)
        {
          dsize = (uint32_t)vector_capacity(_p_vecRate);
        }
        if (iArray) {
                uint32_t* ia = ((uint32_t*)iArray) + *ioffset;
                nrnran123_State** pv = (nrnran123_State**)(&_p_exp_rng);
                nrnran123_getids(*pv, ia, ia+1);

                ia = ia + 2;
                pv = (nrnran123_State**)(&_p_uniform_rng);
                nrnran123_getids(*pv, ia, ia+1);

                ia = ia + 2;
                void* vec = _p_vecRate;
                ia[0] = dsize;

                double *da = dArray + *doffset;
                double *dv;
                if(dsize)
                {
                  dv = vector_vec(vec);
                }
                int iInt;
                for (iInt = 0; iInt < dsize; ++iInt)
                {
                  da[iInt] = dv[iInt];
                }

                vec = _p_vecTbins;
                da = dArray + *doffset + dsize;
                if(dsize)
                {
                  dv = vector_vec(vec);
                }
                for (iInt = 0; iInt < dsize; ++iInt)
                {
                  da[iInt] = dv[iInt];
                }
        }
        *ioffset += 5;
        *doffset += 2*dsize;

}

static void bbcore_read(double* dArray, int* iArray, int* doffset, int* ioffset, _threadargsproto_) {
        assert(!_p_exp_rng);
        assert(!_p_uniform_rng);
        assert(!_p_vecRate);
        assert(!_p_vecTbins);
        uint32_t* ia = ((uint32_t*)iArray) + *ioffset;
        nrnran123_State** pv;
        if (ia[0] != 0 || ia[1] != 0)
        {
          pv = (nrnran123_State**)(&_p_exp_rng);
          *pv = nrnran123_newstream(ia[0], ia[1]);
        }

        ia = ia + 2;
        if (ia[0] != 0 || ia[1] != 0)
        {
          pv = (nrnran123_State**)(&_p_uniform_rng);
          *pv = nrnran123_newstream(ia[0], ia[1]);
        }

        ia = ia + 2;
        int dsize = ia[0];
        *ioffset += 5;

        double *da = dArray + *doffset;
        _p_vecRate = vector_new1(dsize);  /* works for dsize=0 */
        double *dv = vector_vec(_p_vecRate);
        int iInt;
        for (iInt = 0; iInt < dsize; ++iInt)
        {
          dv[iInt] = da[iInt];
        }
        *doffset += dsize;

        da = dArray + *doffset;
        _p_vecTbins = vector_new1(dsize);
        dv = vector_vec(_p_vecTbins);
        for (iInt = 0; iInt < dsize; ++iInt)
        {     
          dv[iInt] = da[iInt];
        }
        *doffset += dsize; 
}
 
static int  setTbins ( _threadargsproto_ ) {
   
/*VERBATIM*/
#if !NRNBBCORE
  void** vv;
  vv = &_p_vecTbins;
  *vv = (void*)0;

  if (ifarg(1)) {
    *vv = vector_arg(1);
  }
#endif
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_setTbins(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 setTbins ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
static int  setRate ( _threadargsproto_ ) {
   
/*VERBATIM*/
#if !NRNBBCORE
  void** vv;
  vv = &_p_vecRate;
  *vv = (void*)0;

  if (ifarg(1)) {
    *vv = vector_arg(1);

    int size = vector_capacity(*vv);
    int i;
    double max=0.0;
    double* px = vector_vec(*vv);
    for (i=0;i<size;i++) {
    	if (px[i]>max) max = px[i];
    }
    rmax = max;
  }
#endif
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_setRate(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 setRate ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
static int  update_time ( _threadargsproto_ ) {
   
/*VERBATIM*/
  void* vv; int i, i_prev, size; double* px;
  i = (int)index;
  i_prev = i;

  if (i >= 0) { // are we disabled?
    vv = _p_vecTbins;
    if (vv) {
      size = vector_capacity(vv);
      px = vector_vec(vv);
      /* advance to current tbins without exceeding array bounds */
      while ((i+1 < size) && (t>=px[i+1])) {
	index += 1.;
	i += 1;
      }
      /* did the index change? */
      if (i!=i_prev) {
        /* advance curRate to next vecRate if possible */
        void *vvRate = _p_vecRate;
        if (vvRate && vector_capacity(vvRate)>i) {
          px = vector_vec(vvRate);
          curRate = px[i];
        }
        else curRate = 1.0;
      }

      /* have we hit last bin? ... disable time advancing leaving curRate as it is*/
      if (i==size)
        index = -1.;

    } else { /* no vecTbins, use some defaults */
      rmax = 1.0;
      curRate = 1.0;
      index = -1.; /* no vecTbins ... disable time advancing & Poisson unit rate. */
    }
  }

  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_update_time(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 update_time ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
void _net_receive2 (_NrnThread * _nt, _Memb_list* _ml, int _iml, int _weight_index, double _lflag);
static void _net_receive (Point_process* _pnt, int _weight_index, double _lflag) 
{   _Memb_list* _ml;  int _iml;

   _NrnThread* _nt;
   int _tid = _pnt->_tid;
   _nt = nrn_threads + _tid;

   _ml = _nt->_ml_list[_pnt->_type];
   _iml = _pnt->_i_instance;
   _net_receive2(_nt, _ml, _iml, _weight_index, _lflag);
}

void _net_receive2 (_NrnThread * _nt, _Memb_list* _ml, int _iml, int _weight_index, double _lflag)
{
  
   assert(0); //BRUNO ADDED THIS (need to solve ref to Point_Process * _pnt first
   Point_process* _pnt = NULL;




   double* _p; Datum* _ppvar;  double v; int _cntml_padded, _cntml_actual; double* _args;
   double *_weights = _nt->_weights;
   ThreadDatum* _thread;

   _thread = (ThreadDatum*)0;
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
 _tsav = t; 
#if !NET_RECEIVE_BUFFERING
  if (_lflag == 1. ) {*(_tqitem) = 0;}
#endif
 {
   if ( _lflag  == 0.0 ) {
     update_time ( _threadargs_ ) ;
     generate_next_event ( _threadargs_ ) ;
     if ( t + event < start + duration ) {
       artcell_net_send ( _tqitem, _weight_index, _pnt, t +  event , 0.0 ) ;
       }
     
/*VERBATIM*/
    
          double u = (double)urand(_threadargs_);
	  /*printf("InhPoisson: spike time at time %g urand=%g curRate=%g, rmax=%g, curRate/rmax=%g \n",t, u, curRate, rmax, curRate/rmax);*/
	  if (u<curRate/rmax) {
 net_event ( _pnt, t ) ;
     
/*VERBATIM*/
    
          }
 }
   } 
#if NET_RECEIVE_BUFFERING
#undef t
#define t _nt->_t
#endif
 }

static void initmodel(_threadargsproto_) {
  int _i; double _save;{
 {
   index = 0. ;
   
/*VERBATIM*/
   if (_p_uniform_rng && _p_exp_rng)
   {
     nrnran123_setseq((nrnran123_State*)_p_uniform_rng, 0, 0);
     nrnran123_setseq((nrnran123_State*)_p_exp_rng, 0, 0);
   }

   void *vvTbins = _p_vecTbins;
   double* px;

   if (vvTbins && vector_capacity(vvTbins)>=1) {
     px = vector_vec(vvTbins);
     start = px[0];
     if (start < 0.0) start=0.0;
   }
   else start = 0.0;

   /* first event is at the start 
   TODO: This should draw from a more appropriate dist
   that has the surrogate process starting a t=-inf
   */
   event = start;

   /* set curRate */
   void *vvRate = _p_vecRate;
   px = vector_vec(vvRate);

   /* set rmax */
   rmax = 0.0;
   int i;
   for (i=0;i<vector_capacity(vvRate);i++) {
      if (px[i]>rmax) rmax = px[i];
   }

   if (vvRate && vector_capacity(vvRate)>0) {
     curRate = px[0];
   }
   else {
      curRate = 1.0;
      rmax = 1.0;
   }

 update_time ( _threadargs_ ) ;
   erand ( _threadargs_ ) ;
   generate_next_event ( _threadargs_ ) ;
   if ( t + event < start + duration ) {
     artcell_net_send ( _tqitem, -1, _nt->_vdata[_ppvar[1*_STRIDE]], t +  event , 0.0 ) ;
     }
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

static double _nrn_current(_threadargsproto_, double _v){double _current=0.;v=_v;{
} return _current;
}

#if defined(ENABLE_CUDA_INTERFACE) && defined(_OPENACC)
  void nrn_state_launcher(_NrnThread*, _Memb_list*, int, int);
  void nrn_jacob_launcher(_NrnThread*, _Memb_list*, int, int);
  void nrn_cur_launcher(_NrnThread*, _Memb_list*, int, int);
#endif


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

#if defined(__cplusplus)
} /* extern "C" */
#endif
