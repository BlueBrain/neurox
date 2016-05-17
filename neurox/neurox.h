#pragma once

#include "hpx/hpx.h"
#include <cassert>

//typedefs
typedef hpx_addr_t hpx_t;
typedef unsigned char byte;

//Memory alignment for hpx_gas_allocs
#define NEUROX_HPX_MEM_ALIGNMENT 0

#define THREAD_ID hpx_thread_get_tls_id()

//just renaming a hpx function for convenience
#define hpx_gas_for hpx_count_range_call //TODO use this, equivalent to a hpx_gas_for (where/how?)

//TODO: using cyclic instead of blocked: this function is not implemented yet
#define hpx_gas_calloc_blocked hpx_gas_calloc_cyclic

#define USE_LCO_FUTURE_ARRAY 0 //TODO: Not working for small node count and high neurons count

#include "neurox/datatypes/InputParams.h"
#include "neurox/datatypes/Mechanism.h"
#include "neurox/datatypes/Compartment.h"
#include "neurox/datatypes/Branch.h"
#include "neurox/datatypes/Neuron.h"
#include "neurox/datatypes/Brain.h"
#include "neurox/input/IDataLoader.h"
#include "neurox/input/CoreNeuronDataLoader.h"


//Global variables (defined on the classes' cpp files)
extern InputParams inputParams;
extern Brain brain;
