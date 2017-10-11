#include "neurox/neurox.h"

#include <algorithm>
#include <cstring>
#include <numeric>
#include <set>

using namespace neurox;
using namespace neurox::tools;
using namespace neurox::interpolators;
using namespace neurox::synchronizers;

const char* BackwardEuler::GetString() {
  return "BackwardEuler";
}

const hpx_action_t BackwardEuler::GetInitAction()
{
    return Finitialize;
}

const hpx_action_t BackwardEuler::GetRunAction()
{
    return RunOnNeuron;
}

const hpx_action_t BackwardEuler::GetRunActionLocality()
{
    return RunOnLocality;
}

int BackwardEuler::GetTotalStepsCount()
{
    return input_params_->tstop_/input_params_->dt_;
}

hpx_action_t BackwardEuler::Finitialize = 0;
int BackwardEuler::Finitialize_handler() {
  NEUROX_MEM_PIN(Branch);
  NEUROX_RECURSIVE_BRANCH_ASYNC_CALL(BackwardEuler::Finitialize);
  Finitialize2(local);
#if !defined(NDEBUG)
// Input::Debugger::StepAfterStepFinitialize(local,
// &nrn_threads[local->nt->id]);
#endif
  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t BackwardEuler::RunOnNeuron = 0;
int BackwardEuler::RunOnNeuron_handler() {
  NEUROX_MEM_PIN(Branch);
  NEUROX_RECURSIVE_BRANCH_ASYNC_CALL(BackwardEuler::RunOnNeuron);
  void * steps_ptr;
  synchronizer_->Run(local, steps_ptr);
  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
  return neurox::wrappers::MemoryUnpin(target);
}


hpx_action_t BackwardEuler::RunOnLocality = 0;
int BackwardEuler::RunOnLocality_handler() {
  NEUROX_MEM_PIN(uint64_t);
  assert(input_params_->allreduce_at_locality_);
  assert(input_params_->synchronizer_ == Synchronizers::kSlidingTimeWindow ||
         input_params_->synchronizer_ == Synchronizers::kAllReduce);

  const int locality_neurons_count =
      AllreduceSynchronizer::AllReducesInfo::AllReduceLocality::
          locality_neurons_->size();
  const hpx_t locality_neurons_lco = hpx_lco_and_new(locality_neurons_count);
  const int comm_step_size = neurox::min_delay_steps_;
  const int reductions_per_comm_step =
      AllreduceSynchronizer::AllReducesInfo::reductions_per_comm_step_;
  const int steps_per_reduction = comm_step_size / reductions_per_comm_step;
  int *steps_ptr;
  const int steps = *steps_ptr;

  for (int s = 0; s < steps; s += comm_step_size) {
    for (int r = 0; r < reductions_per_comm_step; r++) {
      if (s >= comm_step_size)  // first comm-window does not wait
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
                 BackwardEuler::RunOnNeuron, locality_neurons_lco,
                 &steps_per_reduction, sizeof(int));
      hpx_lco_wait_reset(locality_neurons_lco);
    }
  }
  hpx_lco_delete_sync(locality_neurons_lco);
  return neurox::wrappers::MemoryUnpin(target);
}

// fadvance_core.c::nrn_fixed_step_thread
void BackwardEuler::Step(Branch * branch) {
  NrnThread * nt = branch->nt_;
  double &t = nt->_t;
  hpx_t spikes_lco = HPX_NULL;

  synchronizer_->StepBegin(branch);

  // cvodestb.cpp::deliver_net_events()
  // netcvode.cpp::NetCvode::check_thresh(NrnThread*)
  if (branch->soma_) {
    // Soma waits for AIS to have threshold V value updated
    floble_t threshold_v;
    HinesSolver::SynchronizeThresholdV(branch, &threshold_v);
    if (branch->soma_->CheckAPthresholdAndTransmissionFlag(threshold_v))
      spikes_lco = branch->soma_->SendSpikes(nt->_t);
  } else if (branch->thvar_ptr_)
    // Axon Initial Segment send threshold  V to parent
    HinesSolver::SynchronizeThresholdV(branch);

  // netcvode.cpp::NetCvode::deliver_net_events()
  t += .5 * branch->nt_->_dt;
  branch->DeliverEvents(t); // delivers events in the first HALF step
  branch->FixedPlayContinuous();
  SetupTreeMatrix(branch);
  HinesSolver::SolveTreeMatrix(branch);

  // update ions currents based on RHS and dI/dV
  second_order_cur(branch->nt_, input_params_->second_order_);

  ////// fadvance_core.c : update()
  HinesSolver::UpdateVoltagesWithRHS(branch);
  // TODO branch can be placed after the next operation

  // update capacitance currents based on RHS and dI/dV
  branch->CallModFunction(Mechanism::ModFunctions::kCurrentCapacitance);

  ////// fadvance_core.::nrn_fixed_step_lastpart()
  // callModFunction(Mechanism::ModFunction::jacob);
  t += .5 * branch->nt_->_dt;
  // TODO: commenting the call below changes nothing
  //(it changes the variables used in the current function only)
  branch->FixedPlayContinuous();
  branch->CallModFunction(Mechanism::ModFunctions::kState);
  branch->CallModFunction(Mechanism::ModFunctions::kAfterSolve);
  branch->CallModFunction(Mechanism::ModFunctions::kBeforeStep);
  branch->DeliverEvents(t); //delivers events in second HALF-step

  // if we are at the output time instant output to file
  if (fmod(t, input_params_->dt_io_) == 0) {
  }

  synchronizer_->StepEnd(branch, spikes_lco);
}

void BackwardEuler::Finitialize2(Branch * branch) {
  floble_t *v = branch->nt_->_actual_v;
  double t = branch->nt_->_t;

  // set up by finitialize.c:nrn_finitialize(): if (setv)
  assert(input_params_->second_order_ < sizeof(char));
  branch->CallModFunction(Mechanism::ModFunctions::kThreadTableCheck);
  branch->InitVecPlayContinous();
  branch->DeliverEvents(t);

  // set up by finitialize.c:nrn_finitialize(): if (setv)
  for (int n = 0; n < branch->nt_->end; n++) v[n] = input_params_->voltage_;

  // the INITIAL blocks are ordered so that mechanisms that write
  // concentrations are after ions and before mechanisms that read
  // concentrations.
  branch->CallModFunction(Mechanism::ModFunctions::kBeforeInitialize);
  branch->CallModFunction(Mechanism::ModFunctions::kInitialize);
  branch->CallModFunction(Mechanism::ModFunctions::kAfterInitialize);

  branch->DeliverEvents(t);
  SetupTreeMatrix(branch);
  branch->CallModFunction(Mechanism::ModFunctions::kBeforeStep);
  branch->DeliverEvents(t);
}

void BackwardEuler::SetupTreeMatrix(Branch * branch) {
  // treeset_core.c::nrn_rhs: Set up Right-Hand-Side
  // of Matrix-Vector multiplication
  HinesSolver::ResetArray(branch, branch->nt_->_actual_rhs);
  HinesSolver::ResetArray(branch, branch->nt_->_actual_d);

  branch->CallModFunction(Mechanism::ModFunctions::kBeforeBreakpoint);
  branch->CallModFunction(Mechanism::ModFunctions::kCurrent);
  HinesSolver::SetupMatrixRHS(branch);

  // treeset_core.c::nrn_lhs: Set up Left-Hand-Side of Matrix-Vector
  // multiplication. calculate left hand side of
  // cm*dvm/dt = -i(vm) + is(vi) + ai_j*(vi_j - vi)
  // cx*dvx/dt - cm*dvm/dt = -gx*(vx - ex) + i(vm) + ax_j*(vx_j - vx)
  // with a matrix so that the solution is of the form [dvm+dvx,dvx] on the
  // right hand side after solving.
  // branch is a common operation for fixed step, cvode, and daspk methods
  // (TODO: Not for BackwardEuler)
  branch->CallModFunction(Mechanism::ModFunctions::kJacob);

  // finitialize.c:nrn_finitialize()->set_tree_matrix_minimal->nrn_rhs
  // (treeset_core.c)
  // now the cap current can be computed because any change to cm
  // by another model has taken effect.
  branch->CallModFunction(Mechanism::ModFunctions::kJacobCapacitance);

  HinesSolver::SetupMatrixDiagonal(branch);
}


void BackwardEuler::RegisterHpxActions() {
  wrappers::RegisterZeroVarAction(BackwardEuler::Finitialize,
                                  BackwardEuler::Finitialize_handler);
  wrappers::RegisterZeroVarAction(BackwardEuler::RunOnNeuron,
                                  BackwardEuler::RunOnNeuron_handler);
  wrappers::RegisterZeroVarAction(BackwardEuler::RunOnLocality,
                                  BackwardEuler::RunOnLocality_handler);
}
