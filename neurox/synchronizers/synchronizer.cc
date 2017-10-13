#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::synchronizers;

Synchronizer* Synchronizer::New(Synchronizers type) {
  switch (type) {
    case Synchronizers::kDebug:
      return new DebugSynchronizer();
    case Synchronizers::kCoreneuron:
      return new CoreneuronSynchronizer();
    case Synchronizers::kAllReduce:
      return new AllreduceSynchronizer();
    case Synchronizers::kSlidingTimeWindow:
      return new SlidingTimeWindowSynchronizer();
    case Synchronizers::kTimeDependency:
      return new TimeDependencySynchronizer();
    default:
      return nullptr;
  }
  return nullptr;
};

SynchronizerMetadata* SynchronizerMetadata::New(Synchronizers type) {
  switch (type) {
    case Synchronizers::kDebug:
      return new DebugSynchronizer::CommunicationBarrier();
    case Synchronizers::kCoreneuron:
      return new CoreneuronSynchronizer::CommunicationBarrier();
    case Synchronizers::kAllReduce:
      return new AllreduceSynchronizer::AllReducesInfo();
    case Synchronizers::kSlidingTimeWindow:
      return new AllreduceSynchronizer::AllReducesInfo();
    case Synchronizers::kTimeDependency:
      return new TimeDependencySynchronizer::TimeDependencies();
    default:
      return nullptr;
  }
  return nullptr;
}


hpx_action_t Synchronizer::InitLocality = 0;
int Synchronizer::InitLocality_handler(const int *synchronizer_id_ptr, const size_t) {
  NEUROX_MEM_PIN(uint64_t);
  delete neurox::synchronizer_;
  Synchronizers synchronizer_id = *(Synchronizers*)synchronizer_id_ptr;
  neurox::synchronizer_ = Synchronizer::New(synchronizer_id);
  neurox::synchronizer_->Init();;
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::RunLocality = 0;
int Synchronizer::RunLocality_handler(const double * tstop_ptr, const size_t) {
  NEUROX_MEM_PIN(uint64_t);

  const int locality_neurons_count =
      AllreduceSynchronizer::AllReducesInfo::AllReduceLocality::
          locality_neurons_->size();
  const hpx_t locality_neurons_lco = hpx_lco_and_new(locality_neurons_count);
  const int comm_steps = interpolators::BackwardEuler::GetMinSynapticDelaySteps();
  const int reductions_per_comm_step =
      AllreduceSynchronizer::AllReducesInfo::reductions_per_comm_step_;
  const int steps_per_reduction = comm_steps / reductions_per_comm_step;
  const int steps = -1;
  const double tstop = *tstop_ptr;

  for (int s = 0; s < steps; s += comm_steps) {
    for (int r = 0; r < reductions_per_comm_step; r++) {
      if (s >= comm_steps)  // first comm-window does not wait
        hpx_lco_wait_reset(AllreduceSynchronizer::AllReducesInfo::
                               AllReduceLocality::allreduce_future_[r]);
      else
        // fixes crash for Synchronizer::ALL when running two hpx-reduce -based
        // synchronizers in a row
        hpx_lco_reset_sync(AllreduceSynchronizer::AllReducesInfo::
                               AllReduceLocality::allreduce_future_[r]);

      hpx_process_collective_allreduce_join(
          AllreduceSynchronizer::AllReducesInfo::AllReduceLocality::
              allreduce_lco_[r],
          AllreduceSynchronizer::AllReducesInfo::AllReduceLocality::
              allreduce_id_[r],
          NULL, 0);

      for (int i = 0; i < locality_neurons_count; i++)
        hpx_call(AllreduceSynchronizer::AllReducesInfo::AllReduceLocality::
                     locality_neurons_->at(i),
                 Synchronizer::RunNeuron, locality_neurons_lco,
                 &steps_per_reduction, sizeof(int));
      hpx_lco_wait_reset(locality_neurons_lco);
    }
  }
  hpx_lco_delete_sync(locality_neurons_lco);
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::RunNeuron = 0;
int Synchronizer::RunNeuron_handler(const double * tstop_ptr, const size_t) {
  NEUROX_MEM_PIN(Branch);
  NEUROX_RECURSIVE_BRANCH_ASYNC_CALL(Synchronizer::RunNeuron);
  const double tstop= *tstop_ptr;
  local->interpolator_->StepTo(local, tstop);
  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::ClearLocality = 0;
int Synchronizer::ClearLocality_handler() {
  NEUROX_MEM_PIN(uint64_t);
  neurox::synchronizer_->Clear();
  delete neurox::synchronizer_;
  NEUROX_MEM_UNPIN;
}

void Synchronizer::RegisterHpxActions() {
  wrappers::RegisterZeroVarAction(Synchronizer::ClearLocality,
                                  Synchronizer::ClearLocality_handler);
  wrappers::RegisterSingleVarAction<double>(Synchronizer::RunLocality,
                                            Synchronizer::RunLocality_handler);
  wrappers::RegisterSingleVarAction<double>(Synchronizer::RunNeuron,
                                            Synchronizer::RunNeuron_handler);
  wrappers::RegisterSingleVarAction<int>(Synchronizer::InitLocality,
                                         Synchronizer::InitLocality_handler);
}
