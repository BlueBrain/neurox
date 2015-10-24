#pragma once

//Coreneuron core datatypes, tools and methods
#include "coreneuron/coreneuron.h"

//typedefs
typedef double floble_t;    ///> float or double (v, matrix values and mechanisms)
typedef double spike_time_t;///> spikes timing unit
typedef int offset_t;       ///> ushort or uint (p vector, nodes indices)
typedef int neuron_id_t;    ///> neuron gids (gid_t or id_t already used by types.h)
                            // (neuron id is past as -1 on input params parcelgid argument)

//Core datatypes
#include "neurox/hpx.h"
#include "neurox/Event.h"
#include "neurox/VecPlayContinuous.h"
#include "neurox/NetCon.h"
#include "neurox/Mechanism.h"
#include "neurox/Branch.h"
#include "neurox/Neuron.h"

//Fixed-step Backward-Euler solver
#include "neurox/solver/HinesSolver.h"

//Algorithms
#include "neurox/algorithms/Algorithm.h"

//CoreNeuron-based input
#include "neurox/input/Compartment.h"
#include "neurox/input/DataLoader.h"
#include "neurox/input/Debugger.h"

//Tools
#include "neurox/tools/Statistics.h"
#include "neurox/tools/LoadBalancing.h"
#include "neurox/tools/CmdLineParser.h"
#include "neurox/tools/Vectorizer.h"
#ifdef USE_TIM_SPTQ
  #include "neurox/tools/sptq_node.h"
  #include "neurox/tools/sptq_queue.hpp"
  #include "neurox/tools/sptq_queue.ipp"
#endif

//hard-coded mechanism types
#define IClamp 7
#define ProbAMPANMDA_EMS 137
#define ProbGABAAB_EMS 139
#define StochKv 151

//Debugging output
//#define PRINT_TIME_DEPENDENCY

#include <cstring> //std::memset

namespace neurox
{
    extern std::vector<hpx_t> * neurons; ///> hpx address of all neurons

    extern int mechanismsCount; ///> number of mechanisms
    extern neurox::Mechanism ** mechanisms; ///> array to all existing mechanisms
    extern int * mechanismsMap; ///>map of mechanisms offset in 'mechanisms' by 'mechanism type'

    extern tools::CmdLineParser * inputParams; ///> Parameters parsed from command line
    extern algorithms::Algorithm * algorithm; ///> algorithm instance

    extern hpx_action_t Main;            ///> execution starting point (called via hpx_run)
    extern hpx_action_t Clear;           ///> clears all memory utilised including neurons, branches and mechanisms information
    extern hpx_action_t SetMechanisms;   ///> Initializes Mechanisms
    extern hpx_action_t SetMechanismsGlobalVars; ///> sets nrn_ion_global_map and nrn_ion_global_map_size;
    extern hpx_action_t SetAlgorithmVariables;    ///> changes inputParams->algorithm to new value

    Mechanism * GetMechanismFromType(int type); ///> returns mechanisms of type 'type'
    void SetMechanisms2(int count, Mechanism* mechs, int * dependencies, int * successors, char * syms);
    void DebugMessage(const char * str);
    void RunAlgorithm(algorithms::AlgorithmType algorithm);

    static int Main_handler();
    static int Clear_handler();
    static int SetMechanisms_handler (const int nargs, const void *args[], const size_t[]);
    static int SetMechanismsGlobalVars_handler (const int nargs, const void *args[], const size_t[]);
    static int SetAlgorithmVariables_handler(const algorithms::AlgorithmType * algorithm_ptr, const size_t);

    void RegisterHpxActions(); ///> Register all HPX actions
} ;
