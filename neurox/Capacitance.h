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
    void capac_reg_();
    void nrn_cap_jacob(NrnThread* _nt, Memb_list* ml);
    void cap_init(NrnThread* _nt, Memb_list* ml, int type );
    void nrn_capacity_current(NrnThread* _nt, Memb_list* ml);
    void cap_alloc(double*, Datum*, int);
};

};
