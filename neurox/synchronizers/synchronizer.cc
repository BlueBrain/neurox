#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::synchronizers;
using namespace neurox::interpolators;

hpx_t * Synchronizer::locality_neurons_ = nullptr;
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

  //populate local neurons address if not populated
  if (Synchronizer::locality_neurons_==nullptr)
  {
      std::vector<hpx_t> locality_neurons;
      for (int i=0; i<neurox::neurons_count_; i++)
          if (HPX_HERE == HPX_THERE(neurox::neurons_[i]))
              locality_neurons.push_back(neurox::neurons_[i]);
      locality_neurons_ = new hpx_t[locality_neurons.size()];
      std::copy(locality_neurons_, locality_neurons_+locality_neurons.size(), locality_neurons.data());
  }

  //initiate synchronizer
  Synchronizers synchronizer_id = *(Synchronizers*)synchronizer_id_ptr;
  neurox::synchronizer_ = Synchronizer::New(synchronizer_id);
  neurox::synchronizer_->Init();;
  NEUROX_MEM_UNPIN;
}

hpx_action_t Synchronizer::RunLocality = 0;
int Synchronizer::RunLocality_handler(const double * tstop_ptr, const size_t) {
  NEUROX_MEM_PIN(uint64_t);

  assert(locality_neurons_);

  const double tstop = *tstop_ptr - 0.00001;
  const double reduction_interval = 2*input_params_->dt_;

  // fixes crash for Synchronizer::ALL when running two hpx-reduce -based
  // synchronizers in a row
  //locality-reduction-init locality (reset all)
  int r=0;
  hpx_lco_reset_sync(AllreduceSynchronizer::AllReducesInfo::
                         AllReduceLocality::allreduce_future_[r]);

  const hpx_t locality_neurons_lco = hpx_lco_and_new(locality_neurons_count_);

  double step_to_time = -1;
  for (double t=0; t<=tstop; t+= reduction_interval)
  {
      //localit-reduction-begin
      hpx_lco_wait_reset(AllreduceSynchronizer::AllReducesInfo::
                             AllReduceLocality::allreduce_future_[r]);

      hpx_process_collective_allreduce_join(
          AllreduceSynchronizer::AllReducesInfo::AllReduceLocality::
              allreduce_lco_[r],
          AllreduceSynchronizer::AllReducesInfo::AllReduceLocality::
              allreduce_id_[r],
          NULL, 0);


      step_to_time = t+reduction_interval;
      for (int i = 0; i < locality_neurons_count_; i++)
        hpx_call(Synchronizer::locality_neurons_[i],
                 Synchronizer::RunNeuron, locality_neurons_lco,
                 &step_to_time, sizeof(double));

      //locality-reduction- end
      hpx_lco_wait_reset(locality_neurons_lco);
  }


  const int comm_steps = BackwardEuler::GetMinSynapticDelaySteps();
  const int reductions_per_comm_step =
      AllreduceSynchronizer::AllReducesInfo::reductions_per_comm_step_;
  const int steps_per_reduction = comm_steps / reductions_per_comm_step;

  const int total_steps = BackwardEuler::GetTotalSteps();
  for (int s = 0; s < total_steps; s += comm_steps) {
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

      for (int i = 0; i < locality_neurons_count_; i++)
        hpx_call(Synchronizer::locality_neurons_[i],
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
  Interpolator * interpolator = local->interpolator_;
  double & t = local->nt_->_t;
  const double dt_io = input_params_->dt_io_;
  while (t <= tstop)
  {
      synchronizer_->BeforeStep(local);
      hpx_t spikes_lco = interpolator->Step(local);
      synchronizer_->AfterStep(local, spikes_lco);

      if (fmod(t, dt_io) == 0) { //output interval
          //TODO output spikes file
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
