/** When we use solve methods like euler, newton or kinetic schemes,
 *  the state/current updates function need to call solver methods
 *  defined in coreneuron. This is typically done via function pointers.
 *  But for GPU implementation using OpenACC, we can't pass function
 *  pointers.The temporary "workaround" for this was to generate
 *  switch case to select the proper callback function. This is implemented
 *  using python script that look into translated file and generate
 *  _kinderiv.h which has cases for steer functions defined below.
 *  This allows OpenACC to select gpu implementations at compile
 *  time.
 *  \todo: eulerfun/difun are legacy macros and can be replaced with
 *         actual steer function for euler/derivimplicit methods.
 */

#include "coreneuron/mech/cfile/scoplib.h"
#include "coreneuron/mech/mod2c_core_thread.h"
#include "_kinderiv.h"

int derivimplicit_thread(int n, int* slist, int* dlist, DIFUN fun, _threadargsproto_) {
    difun(fun);
    return 0;
}

int nrn_derivimplicit_steer(int fun, _threadargsproto_) {
    switch (fun) { _NRN_DERIVIMPLICIT_CASES }
    return 0;
}

int nrn_euler_steer(int fun, _threadargsproto_) {
    switch (fun) { _NRN_EULER_CASES }
    return 0;
}

int nrn_newton_steer(int fun, _threadargsproto_) {
    switch (fun) { _NRN_DERIVIMPLICIT_NEWTON_CASES }
    return 0;
}

int nrn_kinetic_steer(int fun, SparseObj* so, double* rhs, _threadargsproto_) {
    switch (fun) { _NRN_KINETIC_CASES }
    return 0;
}

// derived from nrn/src/scopmath/euler.c
// updated for aos/soa layout index
int euler_thread(int neqn, int* var, int* der, DIFUN fun, _threadargsproto_) {
    double dt = _nt->_dt;
    int i;

    /* calculate the derivatives */
    eulerfun(fun);

    /* update dependent variables */
    for (i = 0; i < neqn; i++)
        _p[var[i] * _STRIDE] += dt * (_p[der[i] * _STRIDE]);

    return 0;
}
