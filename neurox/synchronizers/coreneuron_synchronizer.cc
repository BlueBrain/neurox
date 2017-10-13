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

void CoreneuronSynchronizer::Launch() {
  assert(0);
  input::Debugger::CompareAllBranches();
}

void CoreneuronSynchronizer::StepBegin(Branch*) { assert(0); }

void CoreneuronSynchronizer::StepEnd(Branch* b, hpx_t) { assert(0); }

void CoreneuronSynchronizer::Run(Branch* b, const void* args) { assert(0); }

hpx_t CoreneuronSynchronizer::SendSpikes(Neuron* neuron, double tt, double) {
  assert(0);
  return HPX_NULL;
}
