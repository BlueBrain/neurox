#pragma once

#include "neurox/Neurox.h"
#include <vector>
#include "coreneuron/nrnoc/membfunc.h" //mod_f_t, mod_alloc_t
#include "coreneuron/nrnoc/nrnoc_ml.h" //ThreadDatum
#include "coreneuron/nrnconf.h" //Datum

namespace Neurox
{
/**
 * @brief The Capacitance namespace
 * Stores unique mechanisms related to functions from the mod files
 */
namespace Capacitance
{
    //from capac.c
    static void cap_alloc(short instancesCount, short dataSize, double * data, short pdataSize, int * pdata, int * nodesIndices);
    static void cap_init(short instancesCount, short dataSize, double * data, short pdataSize, int * pdata, int * nodesIndices);
    static void ion_alloc();
    static void ion_cur(short instancesCount, short dataSize, double * data, short pdataSize, int * pdata, int * nodesIndices);
    static void ion_init(short instancesCount, short dataSize, double * data, short pdataSize, int * pdata, int * nodesIndices);
    static void (*nrn_cap_jacob)(NrnThread*, Memb_list*);
};

};
