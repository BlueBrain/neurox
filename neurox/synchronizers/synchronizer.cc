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

hpx_action_t Synchronizer::CallInitLocality = 0;
int Synchronizer::CallInitLocality_handler(const int* synchronizer_id_ptr,
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

hpx_action_t Synchronizer::CallInitNeuron = 0;
int Synchronizer::CallInitNeuron_handler() {
  NEUROX_MEM_PIN(Branch);
  assert(neurox::synchronizer_);
  neurox::synchronizer_->InitNeuron(local);
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::RunLocality = 0;
int Synchronizer::RunLocality_handler(const double* tstop_ptr, const size_t) {
  NEUROX_MEM_PIN(uint64_t);
  const double reduction_interval =
      synchronizer_->GetLocalityReductionInterval();
  const double tstop = *tstop_ptr;
  double step_to_time = -1;
  for (double t = 0; t <= tstop; t += reduction_interval) {
    synchronizer_->LocalityReduce();
    step_to_time = t + reduction_interval;
    wrappers::CallLocalNeurons(Synchronizer::RunNeuron, &step_to_time, sizeof(double));
  }
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::RunNeuron = 0;
int Synchronizer::RunNeuron_handler(const double* tstop_ptr,
                                    const size_t size) {
  NEUROX_MEM_PIN(Branch);
  NEUROX_RECURSIVE_BRANCH_ASYNC_CALL(Synchronizer::RunNeuron, tstop_ptr, size);
  const double tstop = *tstop_ptr-0.00001;
  Interpolator* interpolator = local->interpolator_;
  const double dt_io = input_params_->dt_io_;
  double tpause = -1;
  hpx_t spikes_lco = HPX_NULL;

  double& t = local->nt_->_t;
  while (t <= tstop) {
    synchronizer_->BeforeSteps(local);
    tpause = synchronizer_->GetMaxStepTime(local);
    spikes_lco = interpolator->StepTo(local, tpause);
    synchronizer_->AfterSteps(local, spikes_lco);
    if (fmod(t, dt_io) == 0) {  // output
    }
  }
  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::CallClearLocality = 0;
int Synchronizer::CallClearLocality_handler() {
  NEUROX_MEM_PIN(uint64_t);
  neurox::synchronizer_->ClearLocality();
  delete neurox::synchronizer_;
  neurox::synchronizer_ = nullptr;
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::CallClearNeuron = 0;
int Synchronizer::CallClearNeuron_handler() {
  NEUROX_MEM_PIN(Branch);
  neurox::synchronizer_->ClearNeuron(local);
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::NeuronInfoDestructor = 0;
int Synchronizer::NeuronInfoDestructor_handler()
{
    NEUROX_MEM_PIN(Branch);
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
  wrappers::RegisterZeroVarAction(Synchronizer::CallClearLocality,
                                  Synchronizer::CallClearLocality_handler);
  wrappers::RegisterZeroVarAction(Synchronizer::CallClearNeuron,
                                  Synchronizer::CallClearNeuron_handler);
  wrappers::RegisterZeroVarAction(Synchronizer::CallInitNeuron,
                                  Synchronizer::CallInitNeuron_handler);
  wrappers::RegisterZeroVarAction(Synchronizer::NeuronInfoDestructor,
                                  Synchronizer::NeuronInfoDestructor_handler);
  wrappers::RegisterSingleVarAction<double>(Synchronizer::RunLocality,
                                            Synchronizer::RunLocality_handler);
  wrappers::RegisterSingleVarAction<double>(Synchronizer::RunNeuron,
                                            Synchronizer::RunNeuron_handler);
  wrappers::RegisterSingleVarAction<int>(
      Synchronizer::CallInitLocality, Synchronizer::CallInitLocality_handler);
  wrappers::RegisterSingleVarAction<int>(Synchronizer::NeuronInfoConstructor,
                                         Synchronizer::NeuronInfoConstructor_handler);
}
