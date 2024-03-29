/*
# =============================================================================
# Copyright (c) 2015 - 2021-2022 Blue Brain Project/EPFL
#
# See top-level LICENSE file for details.
# =============================================================================.
*/
#include "neurox/neurox.h"

#include <algorithm>
#include <cstring>
#include <numeric>
#include <set>

using namespace neurox;
using namespace neurox::tools;
using namespace neurox::interpolators;
using namespace neurox::synchronizers;

const char *BackwardEuler::GetString() { return "BackwardEuler"; }

void BackwardEuler::Init(Branch *b) {
  BackwardEuler::Finitialize2(b);
  b->CallModFunction(Mechanism::ModFunctions::kThreadTableCheck);
#if !defined(NDEBUG)
// Input::Debugger::StepAfterStepFinitialize(local,
// &nrn_threads[local->nt->id]);
#endif
}

void BackwardEuler::StepTo(Branch *b, const double tstop) {
  hpx_t spikes_lco = HPX_NULL;
  while (b->nt_->_t < tstop - 0.000001) {
    // spikes_lco |= BackwardEuler::Step(branch);
    BackwardEuler::Step(b);
    input::Debugger::SingleNeuronStepAndCompare(&nrn_threads[b->nt_->id], b,
                                                input_params_->second_order_);
  }
}

// fadvance_core.c::nrn_fixed_step_thread
hpx_t BackwardEuler::Step(Branch *branch, bool benchmark) {
  NrnThread *nt = branch->nt_;
  double &t = nt->_t;
  hpx_t spikes_lco = HPX_NULL;

  // get previous-step voltages for branching bifurcation
  HinesSolver::UpdateBranchVoltagesWithRHS(branch);

  // cvodestb.cpp::deliver_net_events()
  // netcvode.cpp::NetCvode::check_thresh(NrnThread*)
  if (!benchmark) {
    if (branch->soma_) {
      // perform step-begin operations, eg send updates, wait for dependencies
      synchronizer_->StepSync(branch);

      // Soma waits for AIS to have threshold V value updated
      floble_t threshold_v = HinesSolver::GetAxonInitialSegmentVoltage(branch);

      if (branch->soma_->CheckAPthresholdAndTransmissionFlag(threshold_v))
        spikes_lco = branch->soma_->SendSpikes(nt->_t);
    }
  }

  // netcvode.cpp::NetCvode::deliver_net_events()
  t += .5 * branch->nt_->_dt;
  if (!benchmark) branch->DeliverEvents(t);  // events in 1st half-step
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
  if (!benchmark) branch->DeliverEvents(t);  // events in 2nd half-step

  return spikes_lco;
}

void BackwardEuler::Finitialize2(Branch *branch) {
  floble_t *v = branch->nt_->_actual_v;
  double t = branch->nt_->_t;

  // set up by finitialize.c:nrn_finitialize(): if (setv)
  assert(input_params_->second_order_ < sizeof(char));
  branch->CallModFunction(Mechanism::ModFunctions::kThreadTableCheck);
  branch->InitVecPlayContinous();
  branch->DeliverEvents(t);

  // set up by finitialize.c:nrn_finitialize(): if (setv)
  for (int n = 0; n < branch->nt_->end; n++) v[n] = input_params_->voltage_;
  if (branch->branch_tree_) HinesSolver::InitializeBranchConstants(branch);

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

void BackwardEuler::SetupTreeMatrix(Branch *branch) {
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
