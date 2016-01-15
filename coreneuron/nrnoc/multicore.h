/*
Copyright (c) 2014 EPFL-BBP, All rights reserved.

THIS SOFTWARE IS PROVIDED BY THE BLUE BRAIN PROJECT "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE BLUE BRAIN PROJECT
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef multicore_h
#define multicore_h

#include "coreneuron/nrnconf.h"
#include "coreneuron/nrnoc/membfunc.h"

#if defined(__cplusplus)
class NetCon;
class PreSyn;
extern "C" {
#else
typedef void NetCon;
typedef void PreSyn;
#endif

typedef struct NrnThreadMembList{ /**< patterned after CvMembList in cvodeobj.h */
    struct NrnThreadMembList* next;
    struct Memb_list* ml;
    int index;
    int *dependencies; /**< list of mechanism types that this mechanism depends on*/
    int ndependencies; /**< for scheduling we need to know the dependency count */
} NrnThreadMembList;

typedef struct NrnThreadBAList {
    struct Memb_list* ml; /**< an item in the NrnThreadMembList */
    struct BAMech* bam;
    struct NrnThreadBAList* next;
} NrnThreadBAList;

typedef struct NrnThread {
    double _t;  ///< current simulation time (?)
    double _dt; ///< time step (?)
    double cj; ///<1st or 2nd order solver (?)

	NrnThreadMembList* tml;
	Memb_list** _ml_list;
        Point_process* pntprocs; ///< synapses and artificial cells with and without gid
    PreSyn* presyns; ///< all the output PreSyn with and without gid
        NetCon* netcons;
    double* weights; ///< size n_weight. NetCon.weight_ points into this array.

        int n_pntproc, n_presyn, n_input_presyn, n_netcon, n_weight; ///< only for model_size

        int ncell; /**< analogous to old rootnodecount (? shouldnt this be offset of first node)*/
    int end;    /**< 1 + position of last in v_node array. Now v_node_count. */
    int id; /**< this is nrn_threads[id] */
        int _stop_stepping;
    int n_vecplay; /**< number of instances of VecPlayContinuous */

    size_t _ndata, _nidata, _nvdata; /**< sizes */
    double* _data;   /**< all the other double* and Datum to doubles point into here*/
    int* _idata;     /**< all the Datum to ints index into here */
    void** _vdata;   /**< pointers to Point_Process (?) */
    void** _vecplay; /**< array of instances of VecPlayContinuous */

	double* _actual_rhs;
	double* _actual_d;
	double* _actual_a;
	double* _actual_b;
	double* _actual_v;
	double* _actual_area;
    double* _shadow_rhs; /**< Not pointer into _data. Avoid race for multiple POINT_PROCESS in same compartment */
    double* _shadow_d; /**< Not pointer into _data. Avoid race for multiple POINT_PROCESS in same compartment */
	int* _v_parent_index;
    char* _sp13mat; /**< handle to general sparse matrix */
    struct Memb_list* _ecell_memb_list; /**< normally nil */

    double _ctime; /**< computation time in seconds (using nrnmpi_wtime) */

    NrnThreadBAList* tbl[BEFORE_AFTER_SIZE]; /**< wasteful since almost all empty */

} NrnThread;

extern void nrn_threads_create(int n, int parallel);
extern int nrn_nthread;
extern NrnThread* nrn_threads;
extern void nrn_multithread_job(void*(*)(NrnThread*));
extern void nrn_thread_table_check(void);

extern void nrn_threads_free(void);

#if defined(__cplusplus)
}
#endif

#endif
