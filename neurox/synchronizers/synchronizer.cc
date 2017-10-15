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

hpx_action_t Synchronizer::InitLocalityInfo = 0;
int Synchronizer::InitLocalityInfo_handler(const int* synchronizer_id_ptr,
                                           const size_t) {
  NEUROX_MEM_PIN(uint64_t);

  // delete previous synchronizer (if any)
  delete neurox::synchronizer_;

  // initiate synchronizer
  SynchronizerIds synchronizer_id = *(SynchronizerIds*)synchronizer_id_ptr;
  neurox::synchronizer_ = Synchronizer::New(synchronizer_id);
  neurox::synchronizer_->InitLocality();
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::NeuronInfoConstructor = 0;
int Synchronizer::NeuronInfoConstructor_handler(const int* synchronizer_id_ptr,
                                         const size_t)
{
    NEUROX_MEM_PIN(Branch);
    SynchronizerIds synchronizer_id = *(SynchronizerIds*)synchronizer_id_ptr;
    if (local->soma_)
      local->soma_->synchronizer_neuron_info_ =
          SynchronizerNeuronInfo::New(synchronizer_id);
    NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::InitNeuronInfo = 0;
int Synchronizer::InitNeuronInfo_handler() {
  NEUROX_MEM_PIN(Branch);
  assert(neurox::synchronizer_);
  neurox::synchronizer_->InitNeuron(local);
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::RunLocality = 0;
int Synchronizer::RunLocality_handler(const double* tstop_ptr, const size_t) {
  NEUROX_MEM_PIN(uint64_t);
  const hpx_t locality_neurons_lco =
      hpx_lco_and_new(neurox::locality_neurons_count_);
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

hpx_action_t Synchronizer::RunNeuron = 0;
int Synchronizer::RunNeuron_handler(const double* tstop_ptr,
                                    const size_t size) {
  NEUROX_MEM_PIN(Branch);
  NEUROX_RECURSIVE_BRANCH_ASYNC_CALL(Synchronizer::RunNeuron, tstop_ptr, size);
  const double tstop = *tstop_ptr;
  Interpolator* interpolator = local->interpolator_;

  const double dt_io = input_params_->dt_io_;

  double& t = local->nt_->_t;
  while (t <= tstop) {
    synchronizer_->BeforeSteps(local);
    double tpause = synchronizer_->GetMaxStepTime(local);
    hpx_t spikes_lco = interpolator->StepTo(local, tpause);
    synchronizer_->AfterSteps(local, spikes_lco);
    if (fmod(t, dt_io) == 0) {  // output
    }
  }
  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::ClearLocalityInfo = 0;
int Synchronizer::ClearLocalityInfo_handler() {
  NEUROX_MEM_PIN(uint64_t);
  neurox::synchronizer_->ClearLocality();
  delete neurox::synchronizer_;
  neurox::synchronizer_ = nullptr;
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::ClearNeuronInfo = 0;
int Synchronizer::ClearNeuronInfo_handler() {
  NEUROX_MEM_PIN(Branch);
  neurox::synchronizer_->ClearNeuron(local);
  if (local->soma_) {
    delete local->soma_->synchronizer_neuron_info_;
    local->soma_->synchronizer_neuron_info_ = nullptr;
  }
  NEUROX_MEM_UNPIN;
}

double Synchronizer::GetLocalityReductionInterval()  // virtual
{
  return input_params_->tstop_ + 0.00001;  // i.e. no reduction
}

void Synchronizer::RegisterHpxActions() {
  wrappers::RegisterZeroVarAction(Synchronizer::ClearLocalityInfo,
                                  Synchronizer::ClearLocalityInfo_handler);
  wrappers::RegisterZeroVarAction(Synchronizer::ClearNeuronInfo,
                                  Synchronizer::ClearNeuronInfo_handler);
  wrappers::RegisterZeroVarAction(Synchronizer::InitNeuronInfo,
                                  Synchronizer::InitNeuronInfo_handler);
  wrappers::RegisterSingleVarAction<double>(Synchronizer::RunLocality,
                                            Synchronizer::RunLocality_handler);
  wrappers::RegisterSingleVarAction<double>(Synchronizer::RunNeuron,
                                            Synchronizer::RunNeuron_handler);
  wrappers::RegisterSingleVarAction<int>(
      Synchronizer::InitLocalityInfo, Synchronizer::InitLocalityInfo_handler);
  wrappers::RegisterSingleVarAction<int>(Synchronizer::NeuronInfoConstructor,
                                         Synchronizer::NeuronInfoConstructor_handler);
}
