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

#ifndef nrnoc_ml_h
#define nrnoc_ml_h

#include "coreneuron/nrnconf.h"

typedef union ThreadDatum {
	double val;
	int i;
	double* pval;
	void* _pvoid;
}ThreadDatum;

typedef struct Memb_list {
#if CACHEVEC != 0
    /** nodeindices contains all nodes this extension is responsible for,
     *  ordered according to the matrix. This allows to access the matrix
     *  directly via the nrn_actual_* arrays instead of accessing it in the
     *  order of insertion and via the node-structure, making it more
     *  cache-efficient */
	int *nodeindices;
#endif /* CACHEVEC */
    double* data; ///< offset for the double data in nt._data
    Datum* pdata; ///< offset for the Datum data in nt._nvdata
    ThreadDatum* _thread; /**< thread specific data (when static is no good) */
    int nodecount; /**< actual node count (ie how many nodes have this mechanism ??) */
	int _nodecount_padded;
} Memb_list;

#endif
