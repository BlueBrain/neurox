#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::synchronizers;
using namespace neurox::interpolators;

Synchronizer* Synchronizer::New(SynchronizerIds type) {
  switch (type) {
    case SynchronizerIds::kDebug:
      return new DebugSynchronizer();
    case SynchronizerIds::kCoreneuron:
      return new CoreneuronSynchronizer();
    case SynchronizerIds::kAllReduce:
      return new AllreduceSynchronizer();
    case SynchronizerIds::kSlidingTimeWindow:
      return new SlidingTimeWindowSynchronizer();
    case SynchronizerIds::kTimeDependency:
      return new TimeDependencySynchronizer();
    default:
      return nullptr;
  }
  return nullptr;
};

SynchronizerNeuronInfo* SynchronizerNeuronInfo::New(SynchronizerIds type) {
  switch (type) {
    case SynchronizerIds::kDebug:
      return new DebugSynchronizer::CommunicationBarrier();
    case SynchronizerIds::kCoreneuron:
      return new CoreneuronSynchronizer::CommunicationBarrier();
    case SynchronizerIds::kAllReduce:
      return new AllreduceSynchronizer::AllReduceNeuronInfo(
          AllreduceSynchronizer::kAllReducesCount);
    case SynchronizerIds::kSlidingTimeWindow:
      // same as above, with different number of reduces
      return new AllreduceSynchronizer::AllReduceNeuronInfo(
          SlidingTimeWindowSynchronizer::kAllReducesCount);
    case SynchronizerIds::kTimeDependency:
      return new TimeDependencySynchronizer::TimeDependencies();
    default:
      return nullptr;
  }
  return nullptr;
}

hpx_action_t Synchronizer::InitializeLocality = 0;
int Synchronizer::InitializeLocality_handler(const int* synchronizer_id_ptr,
                                             const size_t) {
  NEUROX_MEM_PIN(uint64_t);

  //delete previous synchronizer (if any)
  delete neurox::synchronizer_;

  // initiate synchronizer
  SynchronizerIds synchronizer_id = *(SynchronizerIds*)synchronizer_id_ptr;
  neurox::synchronizer_ = Synchronizer::New(synchronizer_id);
  neurox::synchronizer_->InitLocality();
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::RunLocality = 0;
int Synchronizer::RunLocality_handler(const double* tstop_ptr, const size_t) {
  NEUROX_MEM_PIN(uint64_t);
  const hpx_t locality_neurons_lco = hpx_lco_and_new(neurox::locality_neurons_count_);
  const double reduction_interval =
      synchronizer_->GetLocalityReductionInterval();
  const double tstop = *tstop_ptr;
  double step_to_time = -1;
  for (double t = 0; t <= tstop; t += reduction_interval) {
    synchronizer_->LocalityReduce();
    step_to_time = t + reduction_interval;
    for (int i = 0; i < neurox::locality_neurons_count_; i++)
      hpx_call(neurox::locality_neurons_[i], Synchronizer::RunNeuron,
               locality_neurons_lco, &step_to_time, sizeof(double));
    hpx_lco_wait_reset(locality_neurons_lco);
  }
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::InitializeNeuron = 0;
int Synchronizer::InitializeNeuron_handler() {
  NEUROX_MEM_PIN(Branch);
  NEUROX_RECURSIVE_BRANCH_ASYNC_CALL(Synchronizer::RunNeuron);
  assert(neurox::synchronizer_);
  neurox::synchronizer_->InitNeuron(local);
  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::RunNeuron = 0;
int Synchronizer::RunNeuron_handler(const double* tstop_ptr, const size_t) {
  NEUROX_MEM_PIN(Branch);
  NEUROX_RECURSIVE_BRANCH_ASYNC_CALL(Synchronizer::RunNeuron);
  const double tstop = *tstop_ptr;
  Interpolator* interpolator = local->interpolator_;

  const double dt_io = input_params_->dt_io_;
  double& t = local->nt_->_t;

  while (t <= tstop) {
    synchronizer_->BeforeStep(local);
    double tpause = synchronizer_->GetMaxStepTime(local);
    hpx_t spikes_lco = interpolator->StepTo(local, tpause);
    synchronizer_->AfterStep(local, spikes_lco);

    if (fmod(t, dt_io) == 0) {  // output
    }
  }
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

double Synchronizer::GetLocalityReductionInterval()  // virtual
{
  return input_params_->tstop_;  // i.e. no reduction
}

void Synchronizer::RegisterHpxActions() {
  wrappers::RegisterZeroVarAction(Synchronizer::ClearLocality,
                                  Synchronizer::ClearLocality_handler);
  wrappers::RegisterZeroVarAction(Synchronizer::InitializeNeuron,
                                  Synchronizer::InitializeNeuron_handler);
  wrappers::RegisterSingleVarAction<double>(Synchronizer::RunLocality,
                                            Synchronizer::RunLocality_handler);
  wrappers::RegisterSingleVarAction<double>(Synchronizer::RunNeuron,
                                            Synchronizer::RunNeuron_handler);
  wrappers::RegisterSingleVarAction<int>(
      Synchronizer::InitializeLocality,
      Synchronizer::InitializeLocality_handler);
}
