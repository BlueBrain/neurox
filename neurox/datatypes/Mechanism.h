#pragma once

#include "neurox/Neurox.h"
#include <vector>
#include "coreneuron/nrnoc/membfunc.h" //mod_f_t, mod_alloc_t
#include "coreneuron/nrnoc/nrnoc_ml.h" //ThreadDatum
#include "coreneuron/nrnconf.h" //Datum


namespace Neurox
{

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

//    typedef void (*mod_alloc_t)(double*, Datum*, int);
//    typedef void (*mod_f_t)(struct NrnThread*, Memb_list*, int);
//    typedef void (*pnt_receive_t)(Point_process*, double*, double);

    //from memb_func.h (before after functions not used on BBP models)
    enum Function
    {
        current=0,
        jacob=1,
        state=2,
        initialize=3,
        before_initialize=4,
        after_initialize=5,
        before_breakpoint=6,
        after_solve=7,
        before_step=8,
        alloc=9,
        functionsCount=10, //number of mod_f_t functions: +1 than previous
        tableCheck=11,
        setData=12,
        capJacob=13,
        destructor=14
    };

    mod_f_t functions[Function::functionsCount]; ///>mechanism functions (with mod_f_t type)

    ///mechanism type -specific functions
    void (*thread_mem_init)(ThreadDatum*); ///> memory allocation
    void (*thread_cleanup)(ThreadDatum*); ///> cleanup
    void (*thread_table_check)(int, int, double*, Datum*, ThreadDatum*, void*, int); ///> Look up tables for mechanisms
    void (*setdata)(double*, Datum*); ///> sets Data??
    void (*pnt_receive_t)(Point_process*, double*, double);

    //Other functions (no mod_f_t type)
    mod_alloc_t alloc;
    Pfri destructor;
    //capacitors only method (capac.c)
    void (*nrn_cap_jacob)(NrnThread*, Memb_list*);


  private:
};

};
