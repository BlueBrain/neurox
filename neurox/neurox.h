#pragma once

//Coreneuron core datatypes
#include "coreneuron/coreneuron.h"

//typedefs
typedef double floble_t;    ///> float or double (v, matrix values and mechanisms)
typedef double spike_time_t;///> spikes timing unit
typedef int offset_t;       ///> ushort or uint (p vector, nodes indices)
typedef int neuron_id_t;    ///> neuron gids (gid_t or id_t already used by types.h)
                            // (neuron id is past as -1 on input params parcelgid argument)

#ifdef  USE_TIMQ
  #include "misc/tim/node.h"
  #include "misc/tim/sptq_queue.hpp"
  #include "misc/tim/sptq_queue.ipp"
#endif

//Core datatypes
#include "neurox/hpx.h"
#include "neurox/Event.h"
#include "neurox/VecPlayContinuous.h"
#include "neurox/NetCon.h"
#include "neurox/Mechanism.h"
#include "neurox/Branch.h"
#include "neurox/Neuron.h"

//Fixed-step Backward-Euler solver
#include "neurox/input/InputParams.h"
#include "neurox/solver/HinesSolver.h"

//CoreNeuron Input
#include "neurox/input/coreneuron/Compartment.h"
#include "neurox/input/coreneuron/DataLoader.h"
#include "neurox/input/coreneuron/Debugger.h"

//Miscellaneous
#include "neurox/tools/Statistics.h"
#include "neurox/tools/LoadBalancing.h"

#define IClamp 7
#define ProbAMPANMDA_EMS 137
#define ProbGABAAB_EMS 139
#define StochKv 151

//Debugging output
//#define PRINT_TIME_DEPENDENCY

namespace neurox
{
    extern std::vector<hpx_t> * neurons; ///> hpx address of all neurons

    extern int mechanismsCount; ///> number of mechanisms
    extern neurox::Mechanism ** mechanisms; ///> array to all existing mechanisms
    extern int * mechanismsMap; ///>map of mechanisms offset in 'mechanisms' by 'mechanism type'

    extern Input::InputParams *  inputParams; ///> Parameters parsed from command line

    extern hpx_action_t main;            ///> execution starting point (called via hpx_run)
    extern hpx_action_t clear;           ///> clears all memory utilised including neurons, branches and mechanisms information
    extern hpx_action_t setMechanisms;   ///> Initializes Mechanisms
    extern hpx_action_t setMechanismsGlobalVars; ///> sets nrn_ion_global_map and nrn_ion_global_map_size;
    extern hpx_action_t setAlgorithmVariables;    ///> changes inputParams->algorithm to new value

    Mechanism * getMechanismFromType(int type); ///> returns mechanisms of type 'type'
    void setMechanisms2(int count, Mechanism* mechs, int * dependencies, int * successors, char * syms);
    void message(const char * str);
    void runAlgorithm(Algorithm algorithm);

    static int main_handler();
    static int clear_handler();
    static int setMechanisms_handler (const int nargs, const void *args[], const size_t[]);
    static int setMechanismsGlobalVars_handler (const int nargs, const void *args[], const size_t[]);
    static int setAlgorithmVariables_handler(const neurox::Algorithm * algorithm_ptr, const size_t);

    void registerHpxActions();           ///> Register all HPX actions
} ;
