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
 
#define nrn_init _nrn_init__BinReport
#define nrn_cur _nrn_cur__BinReport
#define _nrn_current _nrn_current__BinReport
#define nrn_jacob _nrn_jacob__BinReport
#define nrn_state _nrn_state__BinReport
#define initmodel initmodel__BinReport
#define _net_receive _net_receive__BinReport
#define nrn_state_launcher nrn_state_BinReport_launcher
#define nrn_cur_launcher nrn_cur_BinReport_launcher
#define nrn_jacob_launcher nrn_jacob_BinReport_launcher 
#define TimeStatistics TimeStatistics_BinReport 
 
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
#define tstart _p[1*_STRIDE]
#define tstop _p[2*_STRIDE]
#define ISC _p[3*_STRIDE]
#define _v_unused _p[4*_STRIDE]
 
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
 static double _hoc_AddVar();
 static double _hoc_ExtraMapping();
 static double _hoc_TimeStatistics();
 
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
 "AddVar", _hoc_AddVar,
 "ExtraMapping", _hoc_ExtraMapping,
 "TimeStatistics", _hoc_TimeStatistics,
 0, 0
};
 
#endif /*BBCORE*/
#define AddVar AddVar_BinReport
#define ExtraMapping ExtraMapping_BinReport
 inline double AddVar( _threadargsproto_ );
 inline double ExtraMapping( _threadargsproto_ );
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
 "tstart", "ms",
 "tstop", "ms",
 "ISC", "ms",
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
static void _constructor(Prop*);
#endif
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "6.2.0",
"BinReport",
 "Dt",
 "tstart",
 "tstop",
 "ISC",
 0,
 0,
 0,
 "ptr",
 0};
 
static void nrn_alloc(double* _p, Datum* _ppvar, int _type) {
 
#if 0 /*BBCORE*/
 	/*initialize range parameters*/
 	Dt = 0.1;
 	tstart = 0;
 	tstop = 0;
 	ISC = 0;
 
#endif /* BBCORE */
 
#if 0 /*BBCORE*/
if (!nrn_point_prop_) {_constructor(_prop);}
#endif
 
}
 static void _initlists();
 
#define _psize 5
#define _ppsize 3
 void _BinReports_reg() {
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
   hoc_reg_bbcore_read(_mechtype, bbcore_read);
   hoc_reg_bbcore_write(_mechtype, bbcore_write);
  hoc_register_prop_size(_mechtype, _psize, _ppsize);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "bbcorepointer");
 	hoc_register_var(hoc_scdoub, hoc_vdoub, NULL);
 }
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static inline int TimeStatistics(_threadargsproto_);
 } using namespace coreneuron; 
/*VERBATIM*/
#ifndef CORENEURON_BUILD
#ifndef DISABLE_REPORTINGLIB
#include "reportinglib/Records.h"
#include "reportinglib/isc/iscAPI.h"

extern double* hoc_pgetarg(int iarg);
extern double* getarg(int iarg);
extern char* gargstr(int iarg);
extern int hoc_is_str_arg(int iarg);
extern int ifarg(int iarg);
extern double chkarg(int iarg, double low, double high);
//extern double jimsRubbish;

typedef struct {
	char neuronName_[256];
	char rptName_[512];
}Data;

#endif
#endif

 namespace coreneuron { 
double AddVar ( _threadargsproto_ ) {
   double _lAddVar;
 
/*VERBATIM*/
{
#ifndef CORENEURON_BUILD
#ifndef DISABLE_REPORTINGLIB
	Data** tempdata = (Data**)(&(_p_ptr));
	Data* data = *tempdata;
	if(ifarg(1))
	{
		int *ardata =NULL;
		int size = 1;
		ardata = (int *)hoc_Emalloc(sizeof(int)*size);
		hoc_malchk();
		//printf("Mapping info variable Size=%d\n",size);
		ardata[0]=(int) *getarg(2);
		//printf("records_add_var_with_mapping(data->rptName_=%s, data->neuronName_=%s,hoc_pgetarg(1)=%f,size=%d);",data->rptName_,data->neuronName_,hoc_pgetarg(1),size);
		int numberCell;
		sscanf(data->neuronName_,"a%d",&numberCell);

                if ((int)ISC)
                    isc_add_var_with_mapping(data->rptName_, numberCell, hoc_pgetarg(1), ardata);
                else
                    records_add_var_with_mapping(data->rptName_, numberCell, hoc_pgetarg(1), size, ardata);
	}
#endif
#endif
}
 
return _lAddVar;
 }
 
#if 0 /*BBCORE*/
 
static double _hoc_AddVar(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  AddVar ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
double ExtraMapping ( _threadargsproto_ ) {
   double _lExtraMapping;
 
/*VERBATIM*/
{
#ifndef CORENEURON_BUILD
#ifndef DISABLE_REPORTINGLIB
	Data** tempdata = (Data**)(&(_p_ptr));
	Data* data = *tempdata;
	int *ardata =NULL;
	int i=1;
	while(ifarg(i))
	{
		i++;
	}
	int size = i-1;
	//printf("Adding extra mapping. Size=%d\n",size);
	if(size>0)
	{
		ardata = (int *)hoc_Emalloc(sizeof(int)*size);
		hoc_malchk();
		for(i=0;i<size;i++)
		{
			ardata[i]=(int) *getarg(i+1);
			//printf("Adding %d in position %d\n",ardata[i],i);
		}
		int numberCell;
		sscanf(data->neuronName_,"a%d",&numberCell);
		records_extra_mapping(data->rptName_, numberCell,size,ardata);
	}
#endif
#endif
}
 
return _lExtraMapping;
 }
 
#if 0 /*BBCORE*/
 
static double _hoc_ExtraMapping(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  ExtraMapping ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 
static int  TimeStatistics ( _threadargsproto_ ) {
   
/*VERBATIM*/
{
#ifndef CORENEURON_BUILD
#ifndef DISABLE_REPORTINGLIB
         if ((int)ISC == 0)
  	     records_time_data();
#endif
#endif
}
  return 0; }
 
#if 0 /*BBCORE*/
 
static double _hoc_TimeStatistics(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; ThreadDatum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 TimeStatistics ( _threadargs_ );
 return(_r);
}
 
#endif /*BBCORE*/
 } using namespace coreneuron; 
/*VERBATIM*/
/** not executed in coreneuron and hence need empty stubs */
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
#ifndef DISABLE_REPORTINGLIB
	Data** tempdata = (Data**)(&(_p_ptr));//???
	Data* data = 0;
	data = (Data*)hoc_Emalloc(sizeof(Data));
	hoc_malchk();
	if(ifarg(2) && hoc_is_str_arg(2) &&
			ifarg(3) && hoc_is_str_arg(3) &&
			ifarg(4) && hoc_is_str_arg(4) &&
			ifarg(5) && ifarg(6) &&
			ifarg(7) && ifarg(8) && ifarg(9) &&
			ifarg(10) &&
			ifarg(11) && hoc_is_str_arg(11) &&
			ifarg(12) &&
			ifarg(13) && hoc_is_str_arg(13)
	)
	{

		sprintf(data->neuronName_,"%s",gargstr(3));
		sprintf(data->rptName_,"%s/%s",gargstr(4),gargstr(2));
		//printf("%s\n",data->neuronName_);
		//printf("%s\n",data->rptName_);
		unsigned long gid,vgid;
		gid = (unsigned long) *getarg(5);
		vgid = (unsigned long) *getarg(6);
		//printf("Gid %d , Vgid %d\n",gid, vgid);

		tstart = *getarg(7);

		tstop = *getarg(8);

		Dt = *getarg(9);

		//printf("nRpsps received %f, nRpsps recorded %f \n",*getarg(14),reportingSteps);
		//printf("tstart %f , tend %f , dt %f\n",_tstart,_tend,_dt);
		int sizemapping,extramapping;
		sizemapping = (int)*getarg(10);
		extramapping = (int)*getarg(12);
		//printf("sizeMapping %d , extramapping %d\n",sizemapping, extramapping);
		//printf("Hasta aqui\n"); //up to here

		int numberCell;
		sscanf(data->neuronName_,"a%d",&numberCell);

                if (ifarg(14) && hoc_is_str_arg(14))
                {
                    ISC = 1.;
		    isc_add_report(data->rptName_, numberCell, gid, vgid, tstart, tstop, Dt, gargstr(13), gargstr(14));
                }
                else
                {
                    ISC = 0.;
                    records_add_report(data->rptName_,numberCell,gid, vgid, tstart,tstop , Dt,sizemapping,gargstr(11),extramapping,gargstr(13));
                }

		//printf("_________________ %s Creating\n",data->rptName_);

		*tempdata = data; //makes to data available to other procedures through ptr
	}
	else
	{
		//printf("There is an error creating report\n");
		int i = 1;
		while(ifarg(i))
		{
			if(i==1)
				printf("There is an error creating report\n");
			printf("It has arg %d: ");
			if(hoc_is_str_arg(i))
				printf("%s\n",gargstr(i));
			else
				printf("%d\n",*getarg(i));
			i++;
		}

	}
#else
        static warning_shown = 0;
        if (ifarg(2) && hoc_is_str_arg(2))
        {
          if (warning_shown == 0)
          {
            printf("WARNING: BinReports Constructor(): Trying to create and write report %s while the NEURODAMUS_DISABLE_REPORTINGLIB is set to ON, ignoring... \n", gargstr(2));
            warning_shown++;
          }
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
   
/*VERBATIM*/
#ifndef CORENEURON_BUILD
#ifndef DISABLE_REPORTINGLIB
        if ((int)ISC)
            isc_init_connection();
        else
            records_finish_and_share();
#endif
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
    _v = _vec_v[_nd_idx];
    _PRCELLSTATE_V
 v = _v;
 _PRCELLSTATE_V
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
