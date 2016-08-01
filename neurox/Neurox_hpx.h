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

///hpx wrappers for sync call of a function to all children mechanisms
#define neurox_hpx_recursive_mechanism_sync(mechType, Func, ...) \
    short int dependenciesCount = mechanisms[mechType].dependenciesCount; \
    hpx_addr_t lco =  dependenciesCount > 0 ? hpx_lco_and_new(dependenciesCount) : HPX_NULL; \
    for (short int d=0; d<dependenciesCount; d++) \
    { \
        short int dependency = mechanisms[mechType].dependencies[d]; \
        int e = _hpx_call(target, Func, lco, dependency, sizeof(dependency), \
                          __HPX_NARGS(__VA_ARGS__) , ##__VA_ARGS__) ; \
        assert(e==HPX_SUCCESS); \
    } \
    if (dependenciesCount>0) \
    { \
        hpx_lco_wait(lco); \
        hpx_lco_delete(lco, HPX_NULL); \
    }

///hpx wrappers for sync call of a function to all children branches
#define neurox_hpx_recursive_branch_sync(Func, ...) \
    neurox_hpx_recursive_branch_async_call (Func, __HPX_NARGS(__VA_ARGS__) ) \
    neurox_hpx_recursive_branch_async_wait

//auxiliars for neurox_hpx_register_action
#define neurox_hpx_register_action_0(func) \
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE, func, func##_handler);
#define neurox_hpx_register_action_1(func) \
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED, func, func##_handler, \
    HPX_POINTER, HPX_SIZE_T);
#define neurox_hpx_register_action_2(func) \
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED | HPX_VECTORED, \
    func, func##_handler, HPX_INT, HPX_POINTER, HPX_POINTER);
#define neurox_hpx_register_action_3(func) \
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED | HPX_COMPRESSED, \
    func, func##_handler, HPX_POINTER, HPX_SIZE_T);
#define neurox_hpx_register_action_4(func) \
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED | HPX_VECTORED | HPX_COMPRESSED, \
    func, func##_handler, HPX_INT, HPX_POINTER, HPX_POINTER);

/**
 * shortcut for the declaration of hpx functions.
 * Pass function type and name:
 * 0 : no arguments
 * 1 : one argument
 * 2 : more than one arguments
 * 3 : one argument (compressed)
 * 4 : more than one arguments (compressed)
 */
#define neurox_hpx_register_action(funcType, func) \
    neurox_hpx_register_action_##funcType(func)
