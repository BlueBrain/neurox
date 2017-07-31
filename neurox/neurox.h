#pragma once

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
#include "neurox/algorithms/algorithm_metadata.h"

// auxiliary classes defining events, synapses and mechanisms
#include "neurox/event.h"
#include "neurox/mechanism.h"
#include "neurox/net_con.h"
#include "neurox/vec_play_continuous.h"

// morphology classes (branches and soma)
#include "neurox/branch.h"
#include "neurox/neuron.h"

// Fixed-step Backward-Euler solver
#include "neurox/solver/hines_solver.h"

// Tools
#include "neurox/tools/cmd_line_parser.h"
#include "neurox/tools/load_balancing.h"
#include "neurox/tools/statistics.h"
#include "neurox/tools/vectorizer.h"

// Algorithms
#include "neurox/algorithms/algorithm.h"

// CoreNeuron-based input
#include "neurox/input/compartment.h"
#include "neurox/input/data_loader.h"
#include "neurox/input/debugger.h"

// Debug flags
//#define NEUROX_PRINT_TIME_DEPENDENCY
//#define NDEBUG

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

/// Memory pinning from an hpx memorry address to a local pointer
template <typename T>
bool MemoryPin(T *&local, hpx_t &target) {
  target = hpx_thread_current_target();
  return hpx_gas_try_pin(target, (void **)&local);
}

/// Memory UNpinning from an hpx memorry address to a local pointer
template <typename T>
hpx_status_t MemoryUnpin(hpx_t target) {
  hpx_gas_unpin(target);
  return HPX_SUCCESS;
}

/// call action (with arguments) on all localities
template <typename... Args>
hpx_status_t CallAllLocalities(hpx_action_t f, Args... args) {
  return hpx_bcast_rsync(f, args...);
}

/// count the number of arguments
template <typename... ArgTypes>
inline int CountArgs(ArgTypes... args);

template <typename T, typename... ArgTypes>
inline int CountArgs(T t, ArgTypes... args) {
  return 1 + CountArgs(args...);
}
template <>
inline int CountArgs() {
  return 0;
}

/// calls method (with arguments) on all neurons in neurox::neurons
template <typename... Args>
hpx_status_t CallAllNeurons(hpx_action_t f, Args... args) {
  hpx_t lco = hpx_lco_and_new(neurox::neurons_count);
  int e = HPX_SUCCESS;
  int n = neurox::CountArgs(args...);
  for (size_t i = 0; i < neurox::neurons_count; i++)
    e += _hpx_call(neurox::neurons[i], f, lco, n, args...);
  hpx_lco_wait_reset(lco);
  hpx_lco_delete_sync(lco);
  return e;
}
};
