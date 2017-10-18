#include "neurox/synchronizers/debug_synchronizer.h"

using namespace neurox;
using namespace neurox::synchronizers;
using namespace neurox::interpolators;

DebugSynchronizer::DebugSynchronizer() {}

DebugSynchronizer::~DebugSynchronizer() {}

DebugSynchronizer::CommunicationBarrier::CommunicationBarrier() {
  this->all_spikes_lco_ = HPX_NULL;
}

DebugSynchronizer::CommunicationBarrier::~CommunicationBarrier() {
  if (all_spikes_lco_ != HPX_NULL) hpx_lco_delete_sync(all_spikes_lco_);
}

const SynchronizerIds DebugSynchronizer::GetId() {
  return SynchronizerIds::kDebug;
}

const char* DebugSynchronizer::GetString() { return "DebugSynchronizer"; }

void DebugSynchronizer::InitLocality() {}

void DebugSynchronizer::ClearLocality() {}

void DebugSynchronizer::Launch() {
  /*
int comm_step_size = neurox::min_delay_steps_;
int total_steps = 0;

for (int s = 0; s < total_steps; s += comm_step_size) {
#ifdef NEUROX_TIME_STEPPING_VERBOSE
  if (hpx_get_my_rank() == 0)
    DebugMessage(
        std::string("-- t=" + std::to_string(inputParams->dt * s) + " ms\n")
            .c_str());
#endif

  // Reduction at locality not implemented (debugging only)
  wrappers::CallAllNeurons(BackwardEuler::RunOnNeuron, &comm_step_size,
                                   sizeof(int));

#ifndef NDEBUG
  if (neurox::ParallelExecution())  // if parallel execution... spike exchange
    hpx_bcast_rsync(neurox::input::Debugger::NrnSpikeExchange);
#endif
}
input::Debugger::CompareAllBranches();
*/
}

void DebugSynchronizer::BeforeSteps(Branch*) {}

void DebugSynchronizer::AfterSteps(Branch* b, hpx_t) {
  if (b->soma_)  // end of comm-step (steps is the number of steps per commSize)
  {
    CommunicationBarrier* comm_barrier =
        (CommunicationBarrier*)b->soma_->synchronizer_neuron_info_;
    if (comm_barrier->all_spikes_lco_ != HPX_NULL)  // was set/used once
      hpx_lco_wait(comm_barrier->all_spikes_lco_);  // wait if needed
  }
}

double DebugSynchronizer::GetMaxStep(Branch* b) {
  return b->nt_->_dt;  // single step at a time
}

hpx_t DebugSynchronizer::SendSpikes(Neuron* neuron, double tt, double) {
  CommunicationBarrier* comm_barrier =
      (CommunicationBarrier*)neuron->synchronizer_neuron_info_;

  hpx_t& all_spikes_lco = comm_barrier->all_spikes_lco_;

  if (all_spikes_lco == HPX_NULL)  // first use
    all_spikes_lco = hpx_lco_and_new(neuron->synapses_.size());
  else
    hpx_lco_reset_sync(all_spikes_lco);  // reset to use after

  all_spikes_lco = AllreduceSynchronizer::SendSpikes2(neuron, tt);
  return HPX_NULL;
}
