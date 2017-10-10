#include "neurox/synchronizers/debug_synchronizer.h"

using namespace neurox;
using namespace neurox::synchronizers;

DebugSynchronizer::DebugSynchronizer() {}

DebugSynchronizer::~DebugSynchronizer() {}

DebugSynchronizer::CommunicationBarrier::CommunicationBarrier() {
  this->all_spikes_lco_ = HPX_NULL;
  assert(neurox::min_delay_steps_ <= 4);
}

DebugSynchronizer::CommunicationBarrier::~CommunicationBarrier() {
  if (all_spikes_lco_ != HPX_NULL) hpx_lco_delete_sync(all_spikes_lco_);
}

const Synchronizers DebugSynchronizer::GetId() { return Synchronizers::kDebug; }

const char* DebugSynchronizer::GetString() {
  return "BackwardEulerCoreneuronDebug";
}

void DebugSynchronizer::Init() {
  const int allReducesCount = 0;
  Synchronizer::FixedStepMethodsInit();
  hpx_bcast_rsync(
      AllreduceSynchronizer::AllReducesInfo::SetReductionsPerCommStep,
      &allReducesCount, sizeof(int));
}

void DebugSynchronizer::Clear() {}

double DebugSynchronizer::Launch() {
  int comm_step_size = neurox::min_delay_steps_;
  int total_steps = Synchronizer::GetTotalStepsCount();

  hpx_time_t now = hpx_time_now();
  for (int s = 0; s < total_steps; s += comm_step_size) {
#ifdef NEUROX_TIME_STEPPING_VERBOSE
    if (hpx_get_my_rank() == 0)
      DebugMessage(
          std::string("-- t=" + std::to_string(inputParams->dt * s) + " ms\n")
              .c_str());
#endif

    // Reduction at locality not implemented (this is for debugging
    // only)
    neurox::wrappers::CallAllNeurons(Branch::BackwardEuler, &comm_step_size,
                                     sizeof(int));

#ifndef NDEBUG
    if (neurox::ParallelExecution())  // if parallel execution... spike exchange
      hpx_bcast_rsync(neurox::input::Debugger::NrnSpikeExchange);
#endif
  }
  double elapsed_time = hpx_time_elapsed_ms(now) / 1e3;
  input::Debugger::CompareAllBranches();
  return elapsed_time;
}

void DebugSynchronizer::StepBegin(Branch*) {}

void DebugSynchronizer::StepEnd(Branch* b, hpx_t) {
  input::Debugger::SingleNeuronStepAndCompare(&nrn_threads[b->nt_->id], b,
                                              input_params_->second_order_);
}

void DebugSynchronizer::Run(Branch* b, const void* args) {
  int steps = *(int*)args;
  for (int step = 0; step < steps; step++) b->BackwardEulerStep();
  // Input::Coreneuron::Debugger::stepAfterStepBackwardEuler(local,
  // &nrn_threads[this->nt->id], secondorder); //SMP ONLY

  if (b->soma_)  // end of comm-step (steps is the number of steps per commSize)
  {
    CommunicationBarrier* comm_barrier =
        (CommunicationBarrier*)b->soma_->synchronizer_metadata_;
    if (comm_barrier->all_spikes_lco_ != HPX_NULL)  // was set/used once
      hpx_lco_wait(comm_barrier->all_spikes_lco_);  // wait if needed
  }
}

hpx_t DebugSynchronizer::SendSpikes(Neuron* neuron, double tt, double) {
  CommunicationBarrier* comm_barrier =
      (CommunicationBarrier*)neuron->synchronizer_metadata_;
  if (comm_barrier->all_spikes_lco_ == HPX_NULL)  // first use
    comm_barrier->all_spikes_lco_ = hpx_lco_and_new(neuron->synapses_.size());
  else
    hpx_lco_reset_sync(comm_barrier->all_spikes_lco_);  // reset to use after

  for (Neuron::Synapse*& s : neuron->synapses_)
    // deliveryTime (t+delay) is handled on post-syn side (diff value for every
    // NetCon)
    hpx_call(s->branch_addr_, Branch::AddSpikeEvent,
             comm_barrier->all_spikes_lco_, &neuron->gid_, sizeof(neuron_id_t),
             &tt, sizeof(spike_time_t));

  return HPX_NULL;
}
