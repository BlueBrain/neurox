#pragma once

#include "neurox/neurox.h"
#include <vector>
#include "coreneuron/nrnoc/membfunc.h" //mod_f_t, mod_alloc_t
#include "coreneuron/nrnoc/nrnoc_ml.h" //ThreadDatum
#include "coreneuron/nrnconf.h" //Datum

using namespace std;

/**
 * @brief The Mechanisms class
 * Stores unique mechanisms information and dependencies
 */
class Mechanism
{
  public:
    Mechanism(){};
    ~Mechanism();

    Mechanism(const short int type, const short int dataSize, const short int pdataSize,
              const short int dependenciesCount, const char pntMap,
              const char isArtificial, const int * dependencies,
              const char isIon, const double conci, const double conco, const double charge); //for ions only

    short int type, dataSize, pdataSize, dependenciesCount;
    char pntMap, isArtificial;

    static void registerHpxActions(); ///> Register all HPX actions
    int * dependencies; ///> Id of dependencies mechanisms

    //For ionic mechanisms
    int isIon;
    double conci, conco, charge; //from global_conci, global_conco, global_charge variables

    char name[64]; ///> name of the mechanism (variable memb_func[type].sym in CoreNeuron)

    ///mechanism type -specific functions
    void (*thread_mem_init)(ThreadDatum*); ///> memory allocation
    void (*thread_cleanup)(ThreadDatum*); ///> cleanup
    void (*thread_table_check)(int, int, double*, Datum*, ThreadDatum*, void*, int); ///> Look up tables for mechanisms
    void (*setdata)(double*, Datum*); ///> sets Data??

    mod_alloc_t alloc;
    mod_f_t	current;
    mod_f_t	jacob;
    mod_f_t	state;
    mod_f_t	initialize;
    Pfri destructor;

    //capacitors only method (capac.c)
    void (*nrn_cap_jacob)(NrnThread*, Memb_list*){};


    mod_f_t beforeAfterFunctions[BEFORE_AFTER_SIZE]; //Not used in BBP models

  private:
};
