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

#include <math.h>
#include <string.h>
#include "coreneuron/nrnconf.h"
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrnoc/membdef.h"
#include "coreneuron/nrnoc/nrnoc_decl.h"

#if !defined(LAYOUT)
/* 1 means AoS, >1 means AoSoA, <= 0 means SOA */
#define LAYOUT 1
#endif
#if LAYOUT >= 1
#define _STRIDE LAYOUT
#else
#define _STRIDE _cntml_padded + _iml
#endif

extern void hoc_register_prop_size(int, int, int);

#define	nparm 5
static char *mechanism[] = { /*just a template*/
	"0",
	"na_ion",
	"ena", "nao", "nai", 0,
	"ina", "dina_dv_", 0,
	0
};

double nrn_nernst(), nrn_ghk();
static int na_ion, k_ion, ca_ion; /* will get type for these special ions */

int nrn_is_ion(int type) {
	return (memb_func[type].alloc == ion_alloc);
}

static int ion_global_map_size;
static double** ion_global_map;
double** get_ion_global_map() { return ion_global_map;}
#define global_conci(type) ion_global_map[type][0]
#define global_conco(type) ion_global_map[type][1]
#define global_charge(type) ion_global_map[type][2]

double nrn_ion_charge(int type) {
	return global_charge(type);
}

void ion_reg(const char* name, double valence) {
	int i, mechtype;
	char buf[7][50];
	double val;
#define VAL_SENTINAL -10000.

	Sprintf(buf[0], "%s_ion", name);
	Sprintf(buf[1], "e%s", name);
	Sprintf(buf[2], "%si", name);
	Sprintf(buf[3], "%so", name);
	Sprintf(buf[5], "i%s", name);
	Sprintf(buf[6], "di%s_dv_", name);
	for (i=0; i<7; i++) {
		mechanism[i+1] = buf[i];
	}
	mechanism[5] = (char *)0; /* buf[4] not used above */
	mechtype = nrn_get_mechtype(buf[0]);
	if (memb_func[mechtype].alloc != ion_alloc) {
		register_mech((const char**)mechanism, ion_alloc, ion_cur, (mod_f_t)0, (mod_f_t)0, (mod_f_t)ion_init, -1, 1);
		mechtype = nrn_get_mechtype(mechanism[1]);
		_nrn_layout_reg(mechtype, LAYOUT);
		hoc_register_prop_size(mechtype, nparm, 1 );
		nrn_writes_conc(mechtype, 1);
		if (ion_global_map_size <= mechtype) {
			ion_global_map_size = mechtype + 1;
			ion_global_map = (double**)erealloc(ion_global_map,
				sizeof(double*)*ion_global_map_size);
		}
		ion_global_map[mechtype] = (double*)emalloc(3*sizeof(double));
		Sprintf(buf[0], "%si0_%s", name, buf[0]);
		Sprintf(buf[1], "%so0_%s", name, buf[0]);
		if (strcmp("na", name) == 0) {
			na_ion = mechtype;
			global_conci(mechtype) = DEF_nai;
			global_conco(mechtype) = DEF_nao;
			global_charge(mechtype) = 1.;
		}else if (strcmp("k", name) == 0) {
			k_ion = mechtype;
			global_conci(mechtype) = DEF_ki;
			global_conco(mechtype) = DEF_ko;
			global_charge(mechtype) = 1.;
		}else if (strcmp("ca", name) == 0) {
			ca_ion = mechtype;
			global_conci(mechtype) = DEF_cai;
			global_conco(mechtype) = DEF_cao;
			global_charge(mechtype) = 2.;
		}else{
			global_conci(mechtype) = DEF_ioni;
			global_conco(mechtype) = DEF_iono;
			global_charge(mechtype) = VAL_SENTINAL;
		}			
	}
	val = global_charge(mechtype);
	if (valence != VAL_SENTINAL && val != VAL_SENTINAL && valence != val) {
		fprintf(stderr, "%s ion valence defined differently in\n\
two USEION statements (%g and %g)\n",
			buf[0], valence, global_charge(mechtype));
		nrn_exit(1);
	}else if (valence == VAL_SENTINAL && val == VAL_SENTINAL) {
		fprintf(stderr, "%s ion valence must be defined in\n\
the USEION statement of any model using this ion\n", buf[0]);
		nrn_exit(1);
	}else if (valence != VAL_SENTINAL) {
		global_charge(mechtype) = valence;
	}
}

#define FARADAY 96485.309
#define ktf (1000.*8.3134*(celsius + 273.15)/FARADAY)
double nrn_nernst(ci, co, z) double z, ci, co; {
/*printf("nrn_nernst %g %g %g\n", ci, co, z);*/
	if (z == 0) {
		return 0.;
	}
	if (ci <= 0.) {
		return 1e6;
	}else if (co <= 0.) {
		return -1e6;
	}else{
		return ktf/z*log(co/ci);
	}
}

void nrn_wrote_conc(int type, double* p1, int p2, int it, NrnThread* nt) {
    if (it & 04) {
#if LAYOUT <= 0 /* SoA */
        int _iml = 0;
        int _cntml_padded = nt->_ml_list[type]->_nodecount_padded;
#else
        (void)nt;
#endif
        double* pe = p1 - p2*_STRIDE;
        pe[0] = nrn_nernst(pe[1*_STRIDE], pe[2*_STRIDE], nrn_ion_charge(type));
    }
}

static double efun(double x) {
	if (fabs(x) < 1e-4) {
		return 1. - x/2.;
	}else{
		return x/(exp(x) - 1);
	}
}

double nrn_ghk(double v, double ci, double co, double z) {
	double eco, eci, temp;
	temp = z*v/ktf;
	eco = co*efun(temp);
	eci = ci*efun(-temp);
	return (.001)*z*FARADAY*(eci - eco);
}

#if VECTORIZE
#define erev	pd[0*_STRIDE]	/* From Eion */
#define conci	pd[1*_STRIDE]
#define conco	pd[2*_STRIDE]
#define cur	pd[3*_STRIDE]
#define dcurdv	pd[4*_STRIDE]

/*
 handle erev, conci, conc0 "in the right way" according to ion_style
 default. See nrn/lib/help/nrnoc.help.
ion_style("name_ion", [c_style, e_style, einit, eadvance, cinit])

 ica is assigned
 eca is parameter but if conc exists then eca is assigned
 if conc is nrnocCONST then eca calculated on finitialize
 if conc is STATE then eca calculated on fadvance and conc finitialize
 	with global nai0, nao0

 nernst(ci, co, charge) and ghk(v, ci, co, charge) available to hoc
 and models.
*/

#define iontype ppd[0]	/* how _AMBIGUOUS is to be handled */
/*the bitmap is
03	concentration unused, nrnocCONST, DEP, STATE
04	initialize concentrations
030	reversal potential unused, nrnocCONST, DEP, STATE
040	initialize reversal potential
0100	calc reversal during fadvance
0200	ci being written by a model
0400	co being written by a model
*/

#define charge global_charge(type)
#define conci0 global_conci(type)
#define conco0 global_conco(type)

double nrn_nernst_coef(type) int type; {
	/* for computing jacobian element dconc'/dconc */
	return ktf/charge;
}


/* Must be called prior to any channels which update the currents */
void ion_cur(NrnThread* nt, Memb_list* ml, int type) {
	int _cntml_actual = ml->nodecount;
	int _iml;
	double* pd; Datum* ppd;
	(void)nt; /* unused */
/*printf("ion_cur %s\n", memb_func[type].sym->name);*/
#if LAYOUT == 1 /*AoS*/
	for (_iml = 0; _iml < _cntml_actual; ++_iml) {
	  pd = ml->data + _iml*nparm; ppd = ml->pdata + _iml*1;
#endif
#if LAYOUT == 0 /*SoA*/
	int _cntml_padded = ml->_nodecount_padded;
	pd = ml->data; ppd = ml->pdata;
	for (_iml = 0; _iml < _cntml_actual; ++_iml) {
#endif
#if LAYOUT > 1 /*AoSoA*/
#error AoSoA not implemented.
#endif
		dcurdv = 0.;
		cur = 0.;
		if (iontype & 0100) {
			erev = nrn_nernst(conci, conco, charge);
		}
	};
}

/* Must be called prior to other models which possibly also initialize
	concentrations based on their own states
*/
void ion_init(NrnThread* nt, Memb_list* ml, int type) {
	int _cntml_actual = ml->nodecount;
	int _iml;
	double* pd; Datum* ppd;
	(void)nt; /* unused */
/*printf("ion_init %s\n", memb_func[type].sym->name);*/
#if LAYOUT == 1 /*AoS*/
	for (_iml = 0; _iml < _cntml_actual; ++_iml) {
	  pd = ml->data + _iml*nparm; ppd = ml->pdata + _iml*1;
#endif
#if LAYOUT == 0 /*SoA*/
	int _cntml_padded = ml->_nodecount_padded;
	pd = ml->data; ppd = ml->pdata;
	for (_iml = 0; _iml < _cntml_actual; ++_iml) {
#endif
#if LAYOUT > 1 /*AoSoA*/
#error AoSoA not implemented.
#endif
		if (iontype & 04) {
			conci = conci0;
			conco = conco0;
		}
		if (iontype & 040) {
			erev = nrn_nernst(conci, conco, charge);
		}
	}
}

void ion_alloc() {
	assert(0);
}

void second_order_cur(NrnThread* _nt) {
	extern int secondorder;
	NrnThreadMembList* tml;
	Memb_list* ml;
	int _iml, _cntml_actual;
#if LAYOUT == 0
	int _cntml_padded;
#endif
	int* ni;
	double* pd;
	(void)_nt; /* unused */
  if (secondorder == 2) {
	for (tml = _nt->tml; tml; tml = tml->next) if (memb_func[tml->index].alloc == ion_alloc) {
		ml = tml->ml;
		_cntml_actual = ml->nodecount;
		ni = ml->nodeindices;
#if LAYOUT == 1 /*AoS*/
		for (_iml = 0; _iml < _cntml_actual; ++_iml) {
		  pd = ml->data + _iml*nparm;
#endif
#if LAYOUT == 0 /*SoA*/
		_cntml_padded = ml->_nodecount_padded;
		pd = ml->data;
		for (_iml = 0; _iml < _cntml_actual; ++_iml) {
#endif
#if LAYOUT > 1 /*AoSoA*/
#error AoSoA not implemented.
#endif
			cur += dcurdv * ( VEC_RHS(ni[_iml]) );
		}
	}
   }
}

#endif
