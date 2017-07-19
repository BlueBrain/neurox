/**
 * @file NeuroX_hpx.h
 * Includes:
 * - wrappers around HPX functions that are applicable to our use case;
 * - typedefs for HPX related data types
 * - #defines for memory alignment;
 */

#pragma once

#include "hpx/hpx.h"

//typedefs
typedef hpx_addr_t hpx_t;   ///> hpx address (just rephrased with shorter naming)

//Memory alignment for hpx_gas_allocs
#define NEUROX_MEM_ALIGNMENT (2 * sizeof(double))
#define NEUROX_SOA_PADDING 4

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
    return HPX_THREAD_CONTINUE(Var); \
}

///hpx wrappers for async call a function to all children branches (phase 1 - launch threads)
#define neurox_hpx_recursive_branch_async_call(Func, ...) \
    hpx_addr_t lco_branches = local->branchTree && local->branchTree->branchesCount ? hpx_lco_and_new(local->branchTree->branchesCount) : HPX_NULL; \
    if (local->branchTree) \
    for (int c=0; c<local->branchTree->branchesCount; c++) \
       {_hpx_call(local->branchTree->branches[c], Func, lco_branches, __HPX_NARGS(__VA_ARGS__) , ##__VA_ARGS__);} \

///hpx wrappers for async call a function to all children branches (phase 2 - wait for threads)
#define neurox_hpx_recursive_branch_async_wait \
    if (local->branchTree && local->branchTree->branchesCount) \
    { \
        hpx_lco_wait(lco_branches); \
        hpx_lco_delete(lco_branches, HPX_NULL); \
    }

//hpx wrappers to call methods in all neurons
#define neurox_hpx_call_neurons(Func, ...) \
    hpx_par_for_sync( [&] (int i, void*) -> int \
        {  hpx_call_sync(neurox::neurons->at(i), Func, HPX_NULL, 0,  ##__VA_ARGS__); \
        }, 0, neurox::neurons->size(), NULL);

#define neurox_hpx_call_neurons_lco(Func, LCO, ...) \
{ \
    size_t neurons_size = neurons->size(); \
    for (size_t i=0; i< neurons_size; i++) \
        hpx_call(neurox::neurons->at(i), Func, LCO, ##__VA_ARGS__); \
    hpx_lco_wait_reset(LCO); \
}

//Concatenate preprocessor tokens A and B without expanding macro definitions
#define PPCAT_NX(A, B) A ## B
//Concatenate preprocessor tokens A and B after macro-expanding them.
#define PPCAT(A, B) PPCAT_NX(A, B)

//auxiliars for neurox_hpx_register_action (below)
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
#define neurox_hpx_register_action_5(func) \
    HPX_REGISTER_ACTION(HPX_FUNCTION, 0, func, func##_handler);

//main hpx action registration method
#define neurox_zero_var_action 0                 //0: no arguments
#define neurox_single_var_action 1               //1: one argument
#define neurox_several_vars_action 2             //2: more than one argument
#define neurox_single_var_compressed_action 3    //3: one argument (compressed)
#define neurox_several_vars_compressed_action 4  //4: more than one argument (compressed)
#define neurox_reduce_op_action 5                //5: HPX_FUNCTION for reduce operation

#define neurox_hpx_register_action(funcType, func) \
    PPCAT(neurox_hpx_register_action_,funcType)(func)

