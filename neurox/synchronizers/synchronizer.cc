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
                                                const size_t) {
  NEUROX_MEM_PIN(Branch);
  SynchronizerIds sync_id = *(SynchronizerIds*)synchronizer_id_ptr;
  Neuron* soma = local->soma_;
  /* SynchronizerNeuronInfo may have been created before
   * by Neuron::Neuron(...) constructor */
  if (soma && soma->synchronizer_neuron_info_ == nullptr)
    local->soma_->synchronizer_neuron_info_ =
        SynchronizerNeuronInfo::New(sync_id);
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::CallInitNeuron = 0;
int Synchronizer::CallInitNeuron_handler() {
  NEUROX_MEM_PIN(Branch);
  NEUROX_RECURSIVE_BRANCH_ASYNC_CALL(Synchronizer::CallInitNeuron);
  assert(neurox::synchronizer_);
  neurox::synchronizer_->InitNeuron(local);
  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::RunLocality = 0;
int Synchronizer::RunLocality_handler(const double* tstop_ptr, const size_t) {
  NEUROX_MEM_PIN(uint64_t);
  const double reduction_dt = synchronizer_->GetLocalityReductionInterval();
  if (reduction_dt == 0)  // no locality-reduction
    wrappers::CallLocalNeurons(Synchronizer::RunNeuron, tstop_ptr,
                               sizeof(double));
  else {
    const double tstop = *tstop_ptr;
    const double tstart = tstop - input_params_->tstop_;  // for benchmarks
    double step_to_time = -1;
    for (double t = tstart; t <= tstop - 0.00001; t += reduction_dt) {
      synchronizer_->LocalityReduce();
      step_to_time = t + reduction_dt;
      wrappers::CallLocalNeurons(Synchronizer::RunNeuron, &step_to_time,
                                 sizeof(double));
    }
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
  double tpause = -1;
  hpx_t spikes_lco = HPX_NULL;

  double& t = local->nt_->_t;
  while (t < tstop - 0.00001) {
    synchronizer_->BeforeSteps(local);
    tpause = t + synchronizer_->GetMaxStep(local);
    tpause = std::min(tpause, tstop);
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
  NEUROX_RECURSIVE_BRANCH_ASYNC_CALL(Synchronizer::CallClearNeuron);
  neurox::synchronizer_->ClearNeuron(local);
  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::NeuronInfoDestructor = 0;
int Synchronizer::NeuronInfoDestructor_handler() {
  NEUROX_MEM_PIN(Branch);
  if (local->soma_) {
    delete local->soma_->synchronizer_neuron_info_;
    local->soma_->synchronizer_neuron_info_ = nullptr;
  }
  NEUROX_MEM_UNPIN;
}

/// auxiliar method for CallLocalNeurons
hpx_action_t Synchronizer::CallAllNeuronsAux = 0;
int Synchronizer::CallAllNeuronsAux_handler(const int nargs, const void* args[],
                                            const size_t sizes[]) {
  // TODO clean this
  NEUROX_MEM_PIN(uint64_t);
  hpx_action_t f = *(hpx_action_t*)args[0];  // first arg is action id
  if (nargs == 1)
    wrappers::CallLocalNeurons(f);
  else if (nargs == 2)
    wrappers::CallLocalNeurons(f, args[1], sizes[1]);
  else if (nargs == 3)
    wrappers::CallLocalNeurons(f, args[1], sizes[1], args[2], sizes[2]);
  else if (nargs == 4)
    wrappers::CallLocalNeurons(f, args[1], sizes[1], args[2], sizes[2], args[3],
                               sizes[3]);
  else if (nargs == 5)
    wrappers::CallLocalNeurons(f, args[1], sizes[1], args[2], sizes[2], args[3],
                               sizes[3], args[4], sizes[4]);
  else if (nargs == 6)
    wrappers::CallLocalNeurons(f, args[1], sizes[1], args[2], sizes[2], args[3],
                               sizes[3], args[4], sizes[4], args[5], sizes[5]);
  else {
    assert(0);
  }
  NEUROX_MEM_UNPIN;
}

void Synchronizer::RegisterHpxActions() {
  wrappers::RegisterMultipleVarAction(Synchronizer::CallAllNeuronsAux,
                                      Synchronizer::CallAllNeuronsAux_handler);

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
  wrappers::RegisterSingleVarAction<int>(
      Synchronizer::NeuronInfoConstructor,
      Synchronizer::NeuronInfoConstructor_handler);
}
