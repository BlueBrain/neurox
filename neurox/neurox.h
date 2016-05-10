#pragma once

#include "hpx/hpx.h"

namespace neurox;
using namespace neurox;

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

#include "neurox/nrx_setup.h"
#include "neurox/datatypes/GlobalInfo.h"
#include "neurox/datatypes/Neuron.h"
#include "neurox/datatypes/Branch.h"

//Global variables
extern GlobalInfo * globalInfo;
