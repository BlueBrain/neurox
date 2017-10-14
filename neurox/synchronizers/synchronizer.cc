#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::synchronizers;
using namespace neurox::interpolators;

hpx_t* Synchronizer::locality_neurons_ = nullptr;
int Synchronizer::locality_neurons_count_ = -1;

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
      return new AllreduceSynchronizer::AllReduceNeuronInfo(
          AllreduceSynchronizer::kAllReducesCount);
    case Synchronizers::kSlidingTimeWindow:
      // same as above, with different number of reduces
      return new AllreduceSynchronizer::AllReduceNeuronInfo(
          SlidingTimeWindowSynchronizer::kAllReducesCount);
    case Synchronizers::kTimeDependency:
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
  delete neurox::synchronizer_;

  // populate local neurons address if not populated
  if (Synchronizer::locality_neurons_ == nullptr) {
    std::vector<hpx_t> locality_neurons;
    for (int i = 0; i < neurox::neurons_count_; i++)
      if (HPX_HERE == HPX_THERE(neurox::neurons_[i]))
        locality_neurons.push_back(neurox::neurons_[i]);
    locality_neurons_ = new hpx_t[locality_neurons.size()];
    std::copy(locality_neurons_, locality_neurons_ + locality_neurons.size(),
              locality_neurons.data());
  }

  // initiate synchronizer
  Synchronizers synchronizer_id = *(Synchronizers*)synchronizer_id_ptr;
  neurox::synchronizer_ = Synchronizer::New(synchronizer_id);
  neurox::synchronizer_->InitLocality();
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::RunLocality = 0;
int Synchronizer::RunLocality_handler(const double* tstop_ptr, const size_t) {
  NEUROX_MEM_PIN(uint64_t);
  assert(locality_neurons_);
  const hpx_t locality_neurons_lco = hpx_lco_and_new(locality_neurons_count_);
  const double reduction_interval =
      synchronizer_->GetLocalityReductionInterval();
  const double tstop = *tstop_ptr;
  double step_to_time = -1;
  for (double t = 0; t <= tstop; t += reduction_interval) {
    synchronizer_->LocalityReduce();
    step_to_time = t + reduction_interval;
    for (int i = 0; i < locality_neurons_count_; i++)
      hpx_call(Synchronizer::locality_neurons_[i], Synchronizer::RunNeuron,
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
  wrappers::RegisterSingleVarAction<int>(Synchronizer::InitializeLocality,
                                         Synchronizer::InitializeLocality_handler);
}
