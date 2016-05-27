#pragma once

#include "hpx/hpx.h"
#include <cassert>

//typedefs
typedef hpx_addr_t hpx_t;
typedef unsigned char byte;

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

//defines (should be moved to a struct at some point)
#define ALL_NEURONS -1 //TODO make it static

//Memory alignment for hpx_gas_allocs
#define NEUROX_HPX_MEM_ALIGNMENT 0
#define THREAD_ID hpx_thread_get_tls_id()
#define USE_LCO_FUTURE_ARRAY 0 //TODO: Not working for small node count and high neurons count

#include "neurox/datatypes/InputParams.h"
#include "neurox/datatypes/Synapse.h"
#include "neurox/datatypes/Mechanism.h"
#include "neurox/datatypes/Branch.h"
#include "neurox/datatypes/Neuron.h"
#include "neurox/datatypes/Brain.h"
#include "neurox/input/coreneuron/Compartment.h"
#include "neurox/input/coreneuron/DataLoader.h"


#define DataLoader CoreNeuronDataLoader ///> Class responsible for loading data

//Global variables (defined on the classes' cpp files)
extern InputParams * inputParams;
extern Brain * brain;
