#pragma once

//Coreneuron core datatypes
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
#include "neurox/input/InputParams.h"
#include "neurox/solver/HinesSolver.h"

//CoreNeuron Input
#include "neurox/input/coreneuron/Compartment.h"
#include "neurox/input/coreneuron/DataLoader.h"
#include "neurox/input/coreneuron/Debugger.h"

//Miscellaneous
#include "neurox/misc/Statistics.h"

#define IClamp 7
#define ProbAMPANMDA_EMS 137
#define ProbGABAAB_EMS 139

//#define NDEBUG //(if undefined, compares data output with coreneuron and enables assertions)
//#define PRINT_EVENT
//#define PRINT_TIME_DEPENDENCY

///neurox namespace: contains global information that is copied to all localities
namespace neurox
{
    extern std::vector<hpx_t> * neurons; ///> hpx address of all neurons
    extern std::vector<hpx_t> * myNeurons; ///> hpx address of my neurons

    extern int mechanismsCount; ///> number of mechanisms
    extern neurox::Mechanism ** mechanisms; ///> array to all existing mechanisms
    extern int * mechanismsMap; ///>map of mechanisms offset in 'mechanisms' by 'mechanism type'

    extern Input::InputParams *  inputParams; ///> Parameters parsed from command line

    extern hpx_action_t main;            ///> execution starting point (called via hpx_run)
    extern hpx_action_t clear;           ///> clears all memory utilised including neurons, branches and mechanisms information
    extern hpx_action_t setMechanisms;   ///> Initializes Mechanisms
    extern hpx_action_t setMechanismsGlobalVars; ///> sets nrn_ion_global_map and nrn_ion_global_map_size;

    Mechanism * getMechanismFromType(int type); ///> returns mechanisms of type 'type'

    static int main_handler();
    static int clear_handler();
    static int setMechanisms_handler (const int nargs, const void *args[], const size_t[]);
    static int setMechanismsGlobalVars_handler (const int nargs, const void *args[], const size_t[]);

    void setMechanisms2(int count, Mechanism* mechs, int * dependencies, int * successors, char * syms);
    void registerHpxActions();           ///> Register all HPX actions
} ;
