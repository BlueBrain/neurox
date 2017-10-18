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

const SynchronizerIds CoreneuronSynchronizer::GetId() {
  return SynchronizerIds::kCoreneuron;
}

const char* CoreneuronSynchronizer::GetString() {
  return "CoreneuronSynchronizer";
}

void CoreneuronSynchronizer::InitLocality() { assert(0); }

void CoreneuronSynchronizer::ClearLocality() {}

void CoreneuronSynchronizer::NeuronReduceInit(Branch*) { assert(0); }

void CoreneuronSynchronizer::NeuronReduceEnd(Branch* b, hpx_t) { assert(0); }

double CoreneuronSynchronizer::GetReduceInterval()
{ assert(0); }

hpx_t CoreneuronSynchronizer::SendSpikes(Neuron* neuron, double tt, double) {
  assert(0);
  return HPX_NULL;
}
