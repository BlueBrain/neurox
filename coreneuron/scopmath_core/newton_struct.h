#ifndef newton_struct_h
#define newton_struct_h

#include "coreneuron/mech/mod2c_core_thread.h"

#if defined(__cplusplus)
extern "C" {
#endif

/* avoid incessant alloc/free memory */
typedef struct NewtonSpace {
    int n;
    int n_instance;
    double* delta_x;
    double** jacobian;
    int* perm;
    double* high_value;
    double* low_value;
    double* rowmax;
} NewtonSpace;

#pragma acc routine seq
extern int nrn_crout_thread(NewtonSpace* ns, int n, double** a, int* perm, _threadargsproto_);

#pragma acc routine seq
extern void nrn_scopmath_solve_thread(int n,
                                      double** a,
                                      double* value,
                                      int* perm,
                                      double* delta_x,
                                      int* s,
                                      _threadargsproto_);

#pragma acc routine seq
extern int nrn_newton_thread(NewtonSpace* ns,
                             int n,
                             int* s,
                             NEWTFUN pfunc,
                             double* value,
                             _threadargsproto_);

#pragma acc routine seq
extern void nrn_buildjacobian_thread(NewtonSpace* ns,
                                     int n,
                                     int* s,
                                     NEWTFUN pfunc,
                                     double* value,
                                     double** jacobian,
                                     _threadargsproto_);

extern NewtonSpace* nrn_cons_newtonspace(int n, int n_instance);
extern void nrn_destroy_newtonspace(NewtonSpace* ns);

void nrn_newtonspace_copyto_device(NewtonSpace* ns);

#if defined(__cplusplus)
}
#endif

#endif
