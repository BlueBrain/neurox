/**
 * @file Neurox_hpx.h
 * Includes:
 * - wrappers around HPX functions that are applicable to our use case;
 * - typedefs for HPX related data types
 * - #defines for memory alignment;
 */

#pragma once

#include "hpx/hpx.h"

//typedefs
typedef hpx_addr_t hpx_t;

/* TODO
#ifdef DOUBLE_PRECISION
  typedef double floble;
  typedef HPX_DOUBLE HPX_FLOBLE;
#else
  typedef float floble;
  typedef HPX_FLOAT HPX_FLOBLE;
#endif
*/

//Memory alignment for hpx_gas_allocs
#define NEUROX_HPX_MEM_ALIGNMENT 0
#define THREAD_ID hpx_thread_get_tls_id()

///hpx wrappers for the pin operation
#define neurox_hpx_pin(Type) \
    hpx_t target = hpx_thread_current_target(); \
    Type *local = NULL; \
    if (!hpx_gas_try_pin(target, (void**) &local)) \
        return HPX_RESEND;

///hpx wrappers for the unpin operation
#define neurox_hpx_unpin \
{\
    hpx_gas_unpin(target); \
    return HPX_SUCCESS; \
}

///hpx wrappers for the unpin operation and a return value
#define neurox_hpx_unpin_continue(Var) \
{ \
    hpx_gas_unpin(target); \
    HPX_THREAD_CONTINUE(Var); \
}

///hpx wrappers for async call a function to all children branches (phase 1 - launch threads)
#define neurox_hpx_recursive_branch_async_call(Func, ...) \
    hpx_addr_t lco = local->branchesCount ? hpx_lco_and_new(local->branchesCount) : HPX_NULL; \
    for (int c=0; c<local->branchesCount; c++) \
    { \
        int e = _hpx_call(local->branches[c], Func, lco, __HPX_NARGS(__VA_ARGS__) , ##__VA_ARGS__) ; \
        assert(e==HPX_SUCCESS); \
    }

///hpx wrappers for async call a function to all children branches (phase 2 - wait for threads)
#define neurox_hpx_recursive_branch_async_wait \
    if (local->branchesCount>0) \
    { \
        hpx_lco_wait(lco); \
        hpx_lco_delete(lco, HPX_NULL); \
    }

///hpx wrappers for sync call of a function to all children branches
//(hpx_call(local->branches[c], Func, lco, __VA_ARGS__); does not work)
#define neurox_hpx_recursive_branch_sync(Func, ...) \
    neurox_hpx_recursive_branch_async_call (Func, __HPX_NARGS(__VA_ARGS__) ) \
    neurox_hpx_recursive_branch_async_wait
