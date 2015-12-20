#pragma once

// hard-coded mechanism types
#define IClamp 7
#define ProbAMPANMDA_EMS 137
#define ProbGABAAB_EMS 139
#define StochKv 151

// typedefs
typedef double floble_t;  ///> float or double (v, matrix values and mechanisms)
typedef double spike_time_t;  ///> spikes timing unit
typedef int offset_t;         ///> ushort or uint (p vector, nodes indices)
typedef int neuron_id_t;  ///> neuron gid type (gid_t or id_t already defined)

// Coreneuron basic datatypes, input methods, and mechs functions
#include "coreneuron/coreneuron.h"

// hpx macros and typedefs
#include "neurox/hpx.h"

// Definition of meta data specific to a given algorithm
#include "neurox/algorithms/AlgorithmMetaData.h"

// auxiliary classes defining events, synapses and mechanisms
#include "neurox/Event.h"
#include "neurox/VecPlayContinuous.h"
#include "neurox/NetCon.h"
#include "neurox/Mechanism.h"

// morphology classes (branches and soma)
#include "neurox/Branch.h"
#include "neurox/Neuron.h"

// Fixed-step Backward-Euler solver
#include "neurox/solver/HinesSolver.h"

// CoreNeuron-based input
#include "neurox/input/Compartment.h"
#include "neurox/input/DataLoader.h"
#include "neurox/input/Debugger.h"

// Tools
#include "neurox/tools/CmdLineParser.h"
#include "neurox/tools/LoadBalancing.h"
#include "neurox/tools/Statistics.h"
#include "neurox/tools/Vectorizer.h"

// Algorithms
#include "neurox/algorithms/Algorithm.h"

//#define PRINT_TIME_DEPENDENCY

namespace neurox {

///  hpx address of all neurons
extern hpx_t *neurons;

/// length of neurox::neurons
extern int neurons_count;

/// array to all existing mechanisms
extern neurox::Mechanism **mechanisms;

/// length of neuronx::mechanisms
extern int mechanisms_count;

/// map of mechanisms offset in 'mechanisms' by 'mechanism type'
extern int *mechanisms_map;

/// Parameters parsed from command line
extern tools::CmdLineParser *input_params;

/// algorithm instance
extern algorithms::Algorithm *algorithm;

/// returns mechanism of type 'type'
Mechanism *GetMechanismFromType(int type);

/// printf of a given message only on debug mode
void DebugMessage(const char *str);

/// returns true if program launched in more than one locality
bool ParallelExecution();

///  execution starting point (called via hpx_run)
extern hpx_action_t Main;

/// clears all data including neurons, branches and mechanisms information
extern hpx_action_t Clear;

/// handler of HPX-action neurox::Main
static int Main_handler();

/// handler of HPX-action neurox::Clear
static int Clear_handler();

/// HPX-actions registration
void RegisterHpxActions();
};
