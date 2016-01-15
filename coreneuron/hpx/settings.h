#pragma once
  
#include "hpx/hpx.h"

#include "coreneuron/hpx/datatypes.h"


#define THREAD_ID hpx_thread_get_tls_id()

//just renaming a hpx function for convenience
#define hpx_gas_for hpx_count_range_call //TODO use this, equivalent to a hpx_gas_for (where/how?)

//TODO: using cyclic instead of blocked: this function is not implemented yet
#define hpx_gas_calloc_blocked hpx_gas_calloc_cyclic

#define USE_LCO_FUTURE_ARRAY 0 //TODO: Not working for small node count and high neurons count
 

#include "coreneuron/hpx/wrappers.h"

//Global variables
extern GlobalVars * globalVars;
