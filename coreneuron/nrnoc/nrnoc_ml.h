/*
Copyright (c) 2016, Blue Brain Project
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef nrnoc_ml_h
#define nrnoc_ml_h

#define CACHEVEC 2

#define VINDEX	-1
#define CABLESECTION	1
#define MORPHOLOGY	2
#define CAP	3
#define EXTRACELL	5

#define nrnocCONST 1
#define DEP 2
#define STATE 3	/*See init.c and cabvars.h for order of nrnocCONST, DEP, and STATE */

#define BEFORE_INITIAL 0
#define AFTER_INITIAL 1
#define BEFORE_BREAKPOINT 2
#define AFTER_SOLVE 3
#define BEFORE_STEP 4
#define BEFORE_AFTER_SIZE 5 /* 1 more than the previous */

#define VECTORIZE 1

#if defined(__cplusplus)
class NetCon;
class PreSyn;
extern "C" {
#else
typedef void NetCon;
typedef void PreSyn;
#endif

typedef int Datum;
typedef int (*Pfri)();
typedef char Symbol;

typedef struct { const char* name; double* pdoub; } DoubScal;
typedef struct { const char* name; double* pdoub; int index1; } DoubVec;
typedef struct { const char* name; void (*func)(void); } VoidFunc;

struct NrnThread;
typedef struct NrnThread NrnThread;

/* will go away at some point */
typedef struct Point_process {
    int _i_instance;
    short _type;
    short _tid; /* NrnThread id */
} Point_process;

#if PG_ACC_BUGS
typedef struct ThreadDatum {
	int i;
	double* pval;
	void* _pvoid;
}ThreadDatum;
#else
typedef union ThreadDatum {
	double val;
	int i;
	double* pval;
	void* _pvoid;
}ThreadDatum;
#endif

typedef struct NetReceiveBuffer_t {
	int* _pnt_index;
	int* _weight_index;
	double* _nrb_t;
	double* _nrb_flag;
	int _cnt;
	int _size; /* capacity */
	int _pnt_offset;
	int reallocated; /* if buffere resized/reallocated, needs to be copy to gpu */
}NetReceiveBuffer_t;

typedef struct NetSendBuffer_t {
	int* _sendtype; // net_send, net_event, net_move
	int* _vdata_index;
	int* _pnt_index;
	int* _weight_index;
	double* _nsb_t;
	double* _nsb_flag;
	int _cnt;
	int _size; /* capacity */
	int reallocated; /* if buffer resized/reallocated, needs to be copy to cpu */
}NetSendBuffer_t;

///Mechanisms instances handlers
typedef struct Memb_list {
#if CACHEVEC != 0
    /* nodeindices contains all nodes this extension is responsible for,
     * ordered according to the matrix. This allows to access the matrix
     * directly via the nrn_actual_* arrays instead of accessing it in the
     * order of insertion and via the node-structure, making it more
     * cache-efficient */
    int *nodeindices; ///> array of nodes this instance will be applied to
#endif /* CACHEVEC */
    int* _permute;
    double* data; ///> pointer to NrnThread::_data array with start position of this mechanism's data
    Datum* pdata; ///> pointer to NrnThread::pdata array with start position of this mechanism's pointer data
    ThreadDatum* _thread; /* thread specific data (when static is no good) */
    NetReceiveBuffer_t* _net_receive_buffer;
    NetSendBuffer_t* _net_send_buffer;
    int nodecount; ///> number of nodes (instances) of mechanism on the NrnThread
    int _nodecount_padded;
} Memb_list;

//main datatypes of mechanisms functions
typedef void (*mod_alloc_t)(double*, Datum*, int);
typedef void (*mod_f_t)(struct NrnThread*, Memb_list*, int);
typedef void (*pnt_receive_t)(Point_process*, int, double);
typedef void (*bbcore_read_t)(double*, int*, int*, int*, int, int, double*, Datum*, ThreadDatum*, struct NrnThread*, double);

///Mechanisms functions handlers
typedef struct Memb_func {
    mod_alloc_t alloc; //TODO: NOT USED? (see capac.c)
    mod_f_t	current;
    mod_f_t	jacob;
    mod_f_t	state;
    mod_f_t	initialize;
    Pfri	destructor;	/* only for point processes */
    Symbol	*sym;
    int vectorized;
    int thread_size_; /* how many Datum needed in Memb_list if vectorized */
    void (*thread_mem_init_)(ThreadDatum*); /* after Memb_list._thread is allocated */
    void (*thread_cleanup_)(ThreadDatum*); /* before Memb_list._thread is freed */
    void (*thread_table_check_)(int, int, double*, Datum*, ThreadDatum*, void*, int);
    int is_point;
    void (*setdata_)(double*, Datum*);
    int* dparam_semantics; /* for nrncore writing. */
} Memb_func;

///Linked-list of Before-After functions
typedef struct BAMech {
    mod_f_t f;
    int type;
    struct BAMech* next;
} BAMech;

//exposing capacitor functions from nrnoc/capac.c
extern void cap_alloc(double*, Datum*, int type);
extern void cap_cur(NrnThread*, Memb_list*, int);
extern void cap_init(NrnThread*, Memb_list*, int);
extern void cap_jacob(NrnThread*, Memb_list*, int);

//exposing ion functions from nrnoc/eion.c
extern void ion_cur(NrnThread*, Memb_list*, int);
extern void ion_init(NrnThread*, Memb_list*, int);
extern void second_order_cur(NrnThread* _nt);

//exposing all other mechanisms functions from coreneuron/mech/mod_func.c
extern mod_f_t get_init_function(const char * sym);
//extern mod_f_t get_jacob_function(const char * sym);
extern mod_f_t get_cur_function(const char * sym);
extern mod_f_t get_state_function(const char * sym);
extern mod_f_t get_BA_function(const char * sym, int BA_func_id);

//exposing access to ions table (called directly by some mechs)
extern int nrn_ion_global_map_size;
extern double **nrn_ion_global_map;

//exposing temperature variable (used by eion.c)
extern double celsius;

#if defined(__cplusplus)
}
#endif

#endif
