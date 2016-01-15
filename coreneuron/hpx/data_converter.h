#pragma once

#include "coreneuron/hpx/settings.h"

#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrniv/netcvode.h"

void convert_from_coreneuron_to_hpx_datatypes(double tt, NetCvode* ns, NrnThread* nt);
