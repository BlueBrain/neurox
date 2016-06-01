#pragma once

#include "hpx/hpx.h"

//typedefs
typedef hpx_addr_t hpx_t;
typedef unsigned char byte;

namespace Neurox
{
  public:

    extern int neuronsCount; 	///> total neurons count in the system
    extern hpx_t neuronsAddr; 	///> hpx address of the first position of the neurons array

    extern int mechanismsCount;    ///> number of mechanisms
    extern Mechanism * mechanisms; ///> Unique information per mechanism type

    extern Input::InputParams * inputParams; ///> Parameters parsed from command line (TODO should go away at some point)

    inline static hpx_t getNeuronAddr(int i) const {
        return hpx_addr_add(neuronsAddr, sizeof(Neuron)*i, sizeof(Neuron));
    }; ///> returns hpx address for i-th neuron

    static hpx_action_t setInputParams;	///> Initializes InputParams
    static hpx_action_t setMechanisms;	///> Initializes Mechanisms

  private:

    static int init_handler(const int neuronsCount, const hpx_t neuronsAddr,
                            const Mechanism *mechanisms, const size_t mechanismsCount,
                            const int *mechDependencies) ; ///>HPX constructor

    int setInputParams_handler(const Input::InputParams * inputParams, const size_t size); ///handler for setInputParams
    int setMechanisms_handler(const Mechanism * mechanisms, const size_t count); ///handler for setMechanisms
} ;

///hpx wrappers for the pin operation
#define neurox_hpx_pin(Type) \
        hpx_t target = hpx_thread_current_target(); \
        Type *local = NULL; \
        if (!hpx_gas_try_pin(target, (Type**) &local)) \
            return HPX_RESEND;

///hpx wrappers for the unpin operation
#define neurox_hpx_unpin \
        hpx_gas_unpin(target); \
        return HPX_SUCCESS;

#define neurox_hpx_unpin_continue(Var) \
        hpx_gas_unpin(target); \
        HPX_THREAD_CONTINUE(Var);

#define neurox_hpx_recursive_branch_call(Func) \
    hpx_addr_t lco = local->branchesCount ? hpx_lco_and_new(local->branchesCount) : HPX_NULL; \
    for (int c=0; c<local->branchesCount; c++) \
        hpx_call(local->branches[c], Func, lco); \
    if (local->branchesCount>0) \
    { \
        hpx_lco_wait(lco); \
        hpx_lco_delete(lco, HPX_NULL); \
    }

#define neurox_hpx_recursive_branch_call(Func, Var) \
    hpx_addr_t lco = local->branchesCount ? hpx_lco_and_new(local->branchesCount) : HPX_NULL; \
    for (int c=0; c<local->branchesCount; c++) \
        hpx_call(local->branches[c], Func, lco, Var); \
    if (local->branchesCount>0) \
    { \
        hpx_lco_wait(lco); \
        hpx_lco_delete(lco, HPX_NULL); \
    }

//TODO clean up
#define neurox_hpx_recursive_branch_call(Func, Var1, Var2) \
    hpx_addr_t lco = local->branchesCount ? hpx_lco_and_new(local->branchesCount) : HPX_NULL; \
    for (int c=0; c<local->branchesCount; c++) \
        hpx_call(local->branches[c], Func, lco, Var1, Var2); \
    if (local->branchesCount>0) \
    { \
        hpx_lco_wait(lco); \
        hpx_lco_delete(lco, HPX_NULL); \
    }

//defines (should be moved to a struct at some point)
#define ALL_NEURONS -1 //TODO make it static

//Memory alignment for hpx_gas_allocs
#define NEUROX_HPX_MEM_ALIGNMENT 0
#define THREAD_ID hpx_thread_get_tls_id()

#include "neurox/input/InputParams.h"
#include "neurox/datatypes/Synapse.h"
#include "neurox/datatypes/Mechanism.h"
#include "neurox/datatypes/Branch.h"
#include "neurox/datatypes/Neuron.h"

//Data Loaders
#include "neurox/input/coreneuron/Compartment.h"
#include "neurox/input/coreneuron/DataLoader.h"

//Solvers
#include "neurox/solver/BackwardEuler.h"

#define DataLoader CoreNeuronDataLoader ///> Class responsible for loading data


