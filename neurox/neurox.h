#pragma once

#include "hpx/hpx.h"

// typedefs
typedef double floble_t;  ///> float or double (v, matrix values and mechanisms)
typedef double spike_time_t;  ///> spikes timing unit
typedef int offset_t;         ///> ushort or uint (p vector, nodes indices)
typedef int neuron_id_t;   ///> neuron gid type (gid_t or id_t already defined)
typedef hpx_addr_t hpx_t;  ///> hpx address (just rephrased with shorter naming)

// Coreneuron basic datatypes, input methods, and mechs functions
#include "coreneuron/coreneuron.h"

// auxiliary classes defining events, synapses and mechanisms
#include "neurox/event.h"
#include "neurox/mechanism.h"
#include "neurox/netcon.h"
#include "neurox/vecplay_continuous.h"

// morphology classes (branches and soma)
#include "neurox/branch.h"
#include "neurox/neuron.h"

// Tools
#include "neurox/tools/cmd_line_parser.h"
#include "neurox/tools/load_balancing.h"
#include "neurox/tools/statistics.h"
#include "neurox/tools/vectorizer.h"

// Algorithms
#include "neurox/algorithms/algorithm.h"

// Interpolators
#include "neurox/interpolators/interpolator.h"

// Fixed-step Backward-Euler solver
#include "neurox/solver/hines_solver.h"

// CoreNeuron-based input
#include "neurox/input/compartment.h"
#include "neurox/input/data_loader.h"
#include "neurox/input/debugger.h"

// Debug flags
// #define NEUROX_PRINT_TIME_DEPENDENCY
// #define NDEBUG

namespace neurox {

///  hpx address of all neurons
extern hpx_t *neurons_;

/// length of neurox::neurons
extern int neurons_count_;

/// array to all existing mechanisms
extern neurox::Mechanism **mechanisms_;

/// length of neuronx::mechanisms
extern int mechanisms_count_;

/// map of mechanisms offset in 'mechanisms' by 'mechanism type'
extern int *mechanisms_map_;

/// Parameters parsed from command line
extern tools::CmdLineParser *input_params_;

/// neurons synchronization algorithm instance
extern algorithms::Algorithm *algorithm_;

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

// hpx macros and hpx-wrapperss
// TODO can we move this somewhere else?
#include "neurox/macros.h"
#include "neurox/wrappers.h"
