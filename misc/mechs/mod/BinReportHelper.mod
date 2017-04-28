VERBATIM
#ifndef DISABLE_REPORTINGLIB
#include "reportinglib/Records.h"
#include "reportinglib/iscAPI.h"

extern double* hoc_pgetarg(int iarg);
extern double* getarg(int iarg);
extern char* gargstr(int iarg);
extern int hoc_is_str_arg(int iarg);
extern int ifarg(int iarg);
extern double chkarg(int iarg, double low, double high);
extern double* nrn_recalc_ptr(double*);
extern void nrn_register_recalc_ptr_callback(void (*f)(void));

void refreshPointers() { //callback function to update data locations before runtime
	records_refresh_pointers(nrn_recalc_ptr); //tell bin report library to update its pointers using nrn_recalc_ptr function
        isc_refresh_pointers(nrn_recalc_ptr);
}
#endif

ENDVERBATIM

NEURON {
        ARTIFICIAL_CELL BinReportHelper
}


PARAMETER {
	Dt = .1 (ms)
}


INITIAL {
        net_send(0, 1)
}

NET_RECEIVE(w) {	
VERBATIM
	static int step = 0;

#ifndef DISABLE_REPORTINGLIB
	records_rec(step);

        isc_send_data(step);
#endif
	step++;
ENDVERBATIM
	net_send(Dt, 1)
}

CONSTRUCTOR  {
VERBATIM {
#ifndef DISABLE_REPORTINGLIB
	if(ifarg(1))
		Dt = *getarg(1);
	else
		sprintf(stderr,"Invalid agruments for constructor of BinReportHelper\n");

	records_set_atomic_step(Dt);

        isc_set_sim_dt(Dt);

	nrn_register_recalc_ptr_callback( refreshPointers );
#endif
} 
ENDVERBATIM
}


PROCEDURE make_comm() {
VERBATIM 
{
#ifndef DISABLE_REPORTINGLIB
	records_setup_communicator();
#endif
} 
ENDVERBATIM
}


