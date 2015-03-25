/**
 * @file NeuroX_hpx.h
 * Includes:
 * - wrappers around HPX functions that are applicable to our use case;
 * - typedefs for HPX related data types
 * - #defines for memory alignment;
 */

#pragma once

#include "hpx/hpx.h"

// typedefs
typedef hpx_addr_t hpx_t;  ///> hpx address (just rephrased with shorter naming)

// Memory alignment for hpx_gas_allocs and padding (copied from Coreneuron)
#define NEUROX_MEM_ALIGNMENT_ (2 * sizeof(double))
#define NEUROX_SOA_PADDING_ 4

// Threading and locality ids
#define NEUROX_THREAD_ID_ hpx_thread_get_tls_id()
#define NEUROX_MY_RANK_ hpx_get_my_rank()
#define NEUROX_NUM_RANKS_ hpx_get_num_ranks()

/// hpx wrappers for the pin operation
#define NEUROX_MEM_PIN_(Type)                 \
  hpx_t target = hpx_thread_current_target(); \
  Type *local = NULL;                         \
  if (!hpx_gas_try_pin(target, (void **)&local)) return HPX_RESEND;

/// hpx wrappers for the unpin operation
#define NEUROX_MEM_UNPIN_  \
  {                        \
    hpx_gas_unpin(target); \
    return HPX_SUCCESS;    \
  }

/// hpx wrappers for the unpin operation and a return value
#define NEUROX_MEM_UNPIN_CONTINUE_(Var) \
  {                                     \
    hpx_gas_unpin(target);              \
    return HPX_THREAD_CONTINUE(Var);    \
  }

/// hpx wrappers for async call a function to all children branches (phase 1 -
/// launch threads)
#define NEUROX_RECURSIVE_BRANCH_ASYNC_CALL_(Func, ...)              \
  hpx_addr_t lco_branches =                                         \
      local->branchTree && local->branchTree->branchesCount         \
          ? hpx_lco_and_new(local->branchTree->branchesCount)       \
          : HPX_NULL;                                               \
  if (local->branchTree)                                            \
    for (int c = 0; c < local->branchTree->branchesCount; c++) {    \
      _hpx_call(local->branchTree->branches[c], Func, lco_branches, \
                __HPX_NARGS(__VA_ARGS__), ##__VA_ARGS__);           \
    }

/// hpx wrappers for async call a function to all children branches (phase 2 -
/// wait for threads)
#define NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT_                    \
  if (local->branchTree && local->branchTree->branchesCount) { \
    hpx_lco_wait(lco_branches);                                \
    hpx_lco_delete(lco_branches, HPX_NULL);                    \
  }

// hpx wrappers to call methods in all neurons (no args, or only static args)
#define NEUROX_CALL_ALL_NEURONS_(Func, ...)                                  \
  hpx_par_for_sync(                                                          \
      [&](int i, void *) -> int {                                            \
        hpx_call_sync(neurox::neurons[i], Func, HPX_NULL, 0, ##__VA_ARGS__); \
      },                                                                     \
      0, neurox::neurons_count, NULL);

// hpx wrappers to call methods in all neurons (witg args args)
#define NEUROX_CALL_ALL_NEURONS_LCO_(Func, ...)               \
  {                                                           \
    hpx_t LCO = hpx_lco_and_new(neurox::neurons_count);       \
    for (size_t i = 0; i < neurox::neurons_count; i++)        \
      hpx_call(neurox::neurons[i], Func, LCO, ##__VA_ARGS__); \
    hpx_lco_wait_reset(LCO);                                  \
    hpx_lco_delete_sync(LCO);                                 \
  }

// Concatenate preprocessor tokens A and B without expanding macro definitions
#define NEUROX_PPCAT_NX(A, B) A##B
// Concatenate preprocessor tokens A and B after macro-expanding them.
#define NEUROX_PPCAT(A, B) NEUROX_PPCAT_NX(A, B)

// auxiliars for neurox_hpx_register_action (below)
#define NEUROX_REGISTER_ACTION_0(func) \
  HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE, func, func##_handler);
#define NEUROX_REGISTER_ACTION_1(func)                                   \
  HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED, func, func##_handler, \
                      HPX_POINTER, HPX_SIZE_T);
#define NEUROX_REGISTER_ACTION_2(func)                                  \
  HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED | HPX_VECTORED, func, \
                      func##_handler, HPX_INT, HPX_POINTER, HPX_POINTER);
#define NEUROX_REGISTER_ACTION_3(func)                                    \
  HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED | HPX_COMPRESSED, func, \
                      func##_handler, HPX_POINTER, HPX_SIZE_T);
#define NEUROX_REGISTER_ACTION_4(func)                                      \
  HPX_REGISTER_ACTION(HPX_DEFAULT,                                          \
                      HPX_MARSHALLED | HPX_VECTORED | HPX_COMPRESSED, func, \
                      func##_handler, HPX_INT, HPX_POINTER, HPX_POINTER);
#define NEUROX_REGISTER_ACTION_5(func) \
  HPX_REGISTER_ACTION(HPX_FUNCTION, 0, func, func##_handler);

// main hpx action registration method
#define NEUROX_ACTION_ZERO_VAR_ 0                  // no arguments
#define NEUROX_ACTION_SINGLE_VAR_ 1                // one argument
#define NEUROX_ACTION_MULTIPLE_VARS_ 2             // more than one argument
#define NEUROX_ACTION_COMPRESSED_SINGLE_VAR_ 3     // one argument compressed
#define NEUROX_ACTION_COMPRESSED_MULTIPLE_VARS_ 4  // several arguments comp.
#define NEUROX_ACTION_REDUCE_OP_ 5  // HPX_FUNCTION for reduce operation

#define NEUROX_REGISTER_ACTION_(funcType, func) \
  NEUROX_PPCAT(NEUROX_REGISTER_ACTION_, funcType)(func)
