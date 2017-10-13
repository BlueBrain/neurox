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

void BackwardEuler::Init(Branch *b)
{
    BackwardEuler::Finitialize2(b);
    b->CallModFunction(Mechanism::ModFunctions::kThreadTableCheck);
    #if !defined(NDEBUG)
    // Input::Debugger::StepAfterStepFinitialize(local,
    // &nrn_threads[local->nt->id]);
    #endif
}

int BackwardEuler::GetTotalSteps()
{
    return (input_params_->tstop_+0.00001)/input_params_->dt_;
}

int BackwardEuler::GetMinSynapticDelaySteps()
{
    return (neurox::min_synaptic_delay_+0.00001)/input_params_->dt_;
}

void BackwardEuler::FullStep(Branch *branch)
{
    synchronizer_->BeforeStep(branch);
    hpx_t spikes_lco = branch->interpolator_->Step(branch);
    synchronizer_->AfterStep(branch, spikes_lco);
}

// fadvance_core.c::nrn_fixed_step_thread
hpx_t BackwardEuler::Step(Branch* branch)
{
  NrnThread * nt = branch->nt_;
  double &t = nt->_t;
  hpx_t spikes_lco = HPX_NULL;

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

  return spikes_lco;
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
