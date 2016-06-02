#pragma once

#include "neurox/Neurox.h"
#include <vector>
//#include "coreneuron/nrnoc/membfunc.h" //mod_f_t, mod_alloc_t
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

    //from memb_func.h (before after functions not used on BBP models)
    enum Functions
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
        tableCheck=10,
        setData=11,
        capacityCurrent=12, //not in mod files, it's in capac.c
        capJacob=13,
        NetReceive=14,
        destructor=15,
        functionsCount=16, //number of mod_f_t functions: +1 than previous
    };

    typedef void (*modFunction)( short instancesCount, short dataSize, double * data, short pdataSize, int * pdata, int * nodesIndices);
    modFunction functions[Functions::functionsCount]; ///>mechanism functions (with mod_f_t type)

  private:
};

};
