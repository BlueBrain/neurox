#include "neurox/synchronizers/coreneuron_synchronizer.h"

using namespace neurox;
using namespace neurox::synchronizers;

CoreneuronSynchronizer::CoreneuronSynchronizer() {}

CoreneuronSynchronizer::~CoreneuronSynchronizer() {}

CoreneuronSynchronizer::CommunicationBarrier::CommunicationBarrier() {
  assert(0);
}

CoreneuronSynchronizer::CommunicationBarrier::~CommunicationBarrier() {
  assert(0);
}

const Synchronizers CoreneuronSynchronizer::GetId() {
  return Synchronizers::kCoreneuron;
}

const char* CoreneuronSynchronizer::GetString() {
  return "BackwardEulerCoreneuron";
}

void CoreneuronSynchronizer::Init() {
  assert(0);
}

void CoreneuronSynchronizer::Clear() {}

double CoreneuronSynchronizer::Launch() {
  int comm_step_size = neurox::min_delay_steps_;
  hpx_time_t now = hpx_time_now();
  assert(0);
  double elapsed_time = hpx_time_elapsed_ms(now) / 1e3;
  input::Debugger::CompareAllBranches();
  return elapsed_time;
}

void CoreneuronSynchronizer::StepBegin(Branch*) { assert(0); }

void CoreneuronSynchronizer::StepEnd(Branch* b, hpx_t) { assert(0); }

void CoreneuronSynchronizer::Run(Branch* b, const void* args) { assert(0); }

hpx_t CoreneuronSynchronizer::SendSpikes(Neuron* neuron, double tt, double) {
  assert(0);
  return HPX_NULL;
}
