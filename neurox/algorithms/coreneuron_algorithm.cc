#include "neurox/algorithms/coreneuron_algorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

CoreneuronAlgorithm::CoreneuronAlgorithm() {}

CoreneuronAlgorithm::~CoreneuronAlgorithm() {}

CoreneuronAlgorithm::CommunicationBarrier::CommunicationBarrier() { assert(0); }

CoreneuronAlgorithm::CommunicationBarrier::~CommunicationBarrier() {
  assert(0);
}

const SyncAlgorithms CoreneuronAlgorithm::GetId() {
  return SyncAlgorithms::kCoreneuronMPICollective;
}

const char* CoreneuronAlgorithm::GetString() {
  return "BackwardEulerCoreneuron";
}

void CoreneuronAlgorithm::Init() {
  assert(0);
  Algorithm::FixedStepMethodsInit();
}

void CoreneuronAlgorithm::Clear() {}

double CoreneuronAlgorithm::Launch() {
  int comm_step_size = CoreneuronAlgorithm::CommunicationBarrier::kCommStepSize;
  int total_steps = Algorithm::GetTotalStepsCount();
  hpx_time_t now = hpx_time_now();
  assert(0);
  double elapsed_time = hpx_time_elapsed_ms(now) / 1e3;
  input::Debugger::CompareAllBranches();
  return elapsed_time;
}

void CoreneuronAlgorithm::StepBegin(Branch*) { assert(0); }

void CoreneuronAlgorithm::StepEnd(Branch* b, hpx_t) { assert(0); }

void CoreneuronAlgorithm::Run(Branch* b, const void* args) { assert(0); }

hpx_t CoreneuronAlgorithm::SendSpikes(Neuron* neuron, double tt, double) {
  assert(0);
  return HPX_NULL;
}
