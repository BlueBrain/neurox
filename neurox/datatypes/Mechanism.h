#pragma once

#include "neurox/Neurox.h"
#include <vector>
#include "coreneuron/nrnoc/membfunc.h" //mod_f_t, mod_alloc_t
#include "coreneuron/nrnoc/nrnoc_ml.h" //ThreadDatum
#include "coreneuron/nrnconf.h" //Datum
#include "coreneuron/nrnoc/multicore.h" //NrnThread

namespace NeuroX
{
/**
 * @brief The Mechanisms class
 * Stores the unique metadata of each mechanism
 */
class Mechanism
{
  public:
    Mechanism()=delete;
    ~Mechanism();

    Mechanism(const int type, const short int dataSize, const short int pdataSize,
              const char isArtificial, const char pntMap, const char isIon,
              const short int symLengh = 0, const char * sym = nullptr,
              const char isTopMechanism = 0,
              const short int childrenCount = 0, const int * children = nullptr
              );

    int type;
    short int dataSize, pdataSize, vdataSize, childrenCount;
    short int symLength; ///> length of the name of the mechanism;
    char pntMap, isArtificial;
    char isTopMechanism; ///> wether it can be executed directly or requires other to run prior
    char isIon;
    int * children; ///> Id of dependencies mechanisms

    //For ionic mechanisms
    double conci, conco, charge; //from global_conci, global_conco, global_charge variables
    char *sym; ///> name of the mechanism (variable memb_func[type].sym in CoreNeuron)

    //from memb_func.h (before after functions not used on BBP models)
    Memb_func membFunc;
    mod_f_t BAfunctions[BEFORE_AFTER_SIZE]; ///>mechanism functions (with mod_f_t type)
    pnt_receive_t pnt_receive;
    pnt_receive_t pnt_receive_init;
    bbcore_read_t nrn_bbcore_read;

    enum ModFunction
    {
        //BA functions start here (of size BEFORE_AFTER_SIZE)
        before_initialize=0,
        after_initialize=1,
        before_breakpoint=2,
        after_solve=3,
        before_step=4,
        //memb_func functions start here
        alloc=5,
        current=6,
        state=7,
        jacob=8,
        initialize=9,
        destructor=10,
        threadMemInit=11,
        threadCleanup=12,
        threadTableCheck=13,
        setData=14,
        //capacitance functions start here
        currentCapacitance=15, //not in mod files, it's in capac.c
        jacobCapacitance=16,
    };

    ///Call the NetReceive functions on the mod file, for synapses handling
    void callNetReceiveFunction(
            const void * branch, const Spike * spike,
            const char isInitFunction, const double t, const double dt);

    static void registerHpxActions();
    static hpx_action_t callModFunction; ///> calls MOD functions, and BAMembList (nrn_ba)

private:
    void disableMechFunctions(); ///> sets to NULL all function pointers
    void registerIon();  ///> register ions' mechanisms (ion_reg() in eion.c)
    void registerCapacitance();   ///> register mechanism of type "capacitance"
    void registerBeforeAfterFunctions();   ///> register Before-After functions
    void registerModMechanism(); ///> register mechanisms functions (mod_t_f type)
    static int callModFunction_handler(const int nargs, const void *args[], const size_t[]);

};

};
