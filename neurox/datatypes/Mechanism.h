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

    Mechanism(const short type, const short dataSize, const short pdataSize,
              const short dependenciesCount, const char pntMap,
              const char isArtificial, const int * dependencies,
              const char isIon, const double conci, const double conco, const double charge); //for ions only

    short int type, dataSize, pdataSize, dependenciesCount;
    char pntMap, isArtificial;
    int * dependencies; ///> Id of dependencies mechanisms

    //For ionic mechanisms
    int isIon;
    double conci, conco, charge; //from global_conci, global_conco, global_charge variables

    char name[64]; ///> name of the mechanism (variable memb_func[type].sym in CoreNeuron)

    //from memb_func.h (before after functions not used on BBP models)
    Memb_func membFunc;
    mod_f_t BAfunctions[BEFORE_AFTER_SIZE]; ///>mechanism functions (with mod_f_t type)
    pnt_receive_t pnt_receive;
    pnt_receive_t pnt_receive_init;

    enum ModFunction
    {
        //BA functions start here (of size BEFORE_AFTER_SIZE)
        before_initialize=0,
        after_initialize=1,
        before_breakpoint=2,
        after_solve=3,
        before_step=4,
        //memb_func functions start here
        current=5,
        jacob=6,
        state=7,
        initialize=8,
        alloc=9,
        threadMemInit=10,
        threadCleanup=11,
        threadTableCheck=12,
        setData=13,
        destructor=14,
        //capacitance functions start here
        capacityCurrent=15, //not in mod files, it's in capac.c
        capJacob=16,
    };

    static void registerHpxActions(); ///> Register all HPX actions
    static hpx_action_t setMechanisms; ///> Sets Mechanism

private:
    static int setMechanisms_handler(const Mechanism * mechanisms, const size_t);
};

};
