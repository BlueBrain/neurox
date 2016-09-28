#include "coreneuron/mech/cfile/scoplib.h"
#include "coreneuron/mech/mod2c_core_thread.h"
#include "coreneuron/scopmath_core/errcodes.h"
#define s_(arg) _p[s[arg] * _STRIDE]

#pragma acc routine seq
static int check_state(int, int*, _threadargsproto_);

int _ss_sparse_thread(SparseObj* v,
                      int n,
                      int* s,
                      int* d,
                      double* t,
                      double dt,
                      SPFUN fun,
                      int linflag,
                      _threadargsproto_) {
    int err, i;
    double ss_dt;

    ss_dt = 1e9;
    _modl_set_dt_thread(ss_dt, _nt);

    if (linflag) { /*iterate linear solution*/
        err = sparse_thread(v, n, s, d, t, ss_dt, fun, 0, _threadargs_);
    } else {
#define NIT 7
        i = NIT;
        err = 0;
        while (i) {
            err = sparse_thread(v, n, s, d, t, ss_dt, fun, 1, _threadargs_);
            if (!err) {
                if (check_state(n, s, _threadargs_)) {
                    err = sparse_thread(v, n, s, d, t, ss_dt, fun, 0, _threadargs_);
                }
            }
            --i;
            if (!err) {
                i = 0;
            }
        }
    }

    _modl_set_dt_thread(dt, _nt);
    return err;
}

int _ss_derivimplicit_thread(int n, int* slist, int* dlist, DIFUN fun, _threadargsproto_) {
    int err, i;
    double dtsav;

    dtsav = _modl_get_dt_thread(_nt);
    _modl_set_dt_thread(1e-9, _nt);

    err = derivimplicit_thread(n, slist, dlist, fun, _threadargs_);

    _modl_set_dt_thread(dtsav, _nt);
    return err;
}

static int check_state(int n, int* s, _threadargsproto_) {
    int i, flag;

    flag = 1;
    for (i = 0; i < n; i++) {
        if (s_(i) < -1e-6) {
            s_(i) = 0.;
            flag = 0;
        }
    }
    return flag;
}

void _modl_set_dt_thread(double dt, NrnThread* nt) {
    nt->_dt = dt;
}
double _modl_get_dt_thread(NrnThread* nt) {
    return nt->_dt;
}
