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
    Mechanism();
    ~Mechanism();

    Mechanism(const short int type, const short int dataSize, const short int pdataSize,
              const char isArtificial, const char pntMap, const char isIon,
              const short int symLengh = 0, const char * sym = nullptr,
              const short int dependenciesCount = 0, const short int * dependencies = nullptr);

    short int type, dataSize, pdataSize, dependenciesCount;
    char pntMap, isArtificial;
    short int symLength; ///> length of the name of the mechanism;
    int isIon;
    short int * dependencies; ///> Id of dependencies mechanisms

    //For ionic mechanisms
    double conci, conco, charge; //from global_conci, global_conco, global_charge variables
    char *sym; ///> name of the mechanism (variable memb_func[type].sym in CoreNeuron)

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

    static void registerHpxActions();
    static hpx_action_t callModFunction; ///> calls MOD functions, and BAMembList (nrn_ba)
    static hpx_action_t callNetReceiveFunction; ///> calls NetReceive Functions

private:
    void registerIonicCharges();  ///> register ionic charges (ion_reg() in eion.c)
    void registerBAFunctions();   ///> register Before-After functions
    void registerMechFunctions(); ///> register mechanisms functions (mod_t_f type)
    static int callModFunction_handler(const int nargs, const void *args[], const size_t[]);
    static int callNetReceiveFunction_handler(const int nargs, const void *args[], const size_t[]);
};

};
