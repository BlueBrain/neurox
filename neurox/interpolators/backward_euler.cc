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

int BackwardEuler::GetTotalSteps() {
  return (input_params_->tstop_ + 0.00001) / input_params_->dt_;
}

int BackwardEuler::GetMinSynapticDelaySteps() {
  return (neurox::min_synaptic_delay_ + 0.00001) / input_params_->dt_;
}

hpx_t BackwardEuler::StepTo(Branch *b, const double tstop) {
  hpx_t spikes_lco = HPX_NULL;
  while (b->nt_->_t < tstop - 0.000001) {
    // spikes_lco |= BackwardEuler::Step(b);
    // hpx_t new_spikes_lco = BackwardEuler::Step(b);
    hpx_t new_spikes_lco = BackwardEuler::StepParallel(b);
    input::Debugger::SingleNeuronStepAndCompare(&nrn_threads[b->nt_->id], b,
                                                input_params_->second_order_);
    if (new_spikes_lco) {
      // make sure only one AP occurred in between
      assert(!spikes_lco);
      spikes_lco = new_spikes_lco;
    }
  }
  return spikes_lco;
}

// fadvance_core.c::nrn_fixed_step_thread
hpx_t BackwardEuler::Step(Branch *branch, bool benchmark) {
  NrnThread *nt = branch->nt_;
  double &t = nt->_t;
  const double dt = nt->_dt;
  hpx_t spikes_lco = HPX_NULL;

  // cvodestb.cpp::deliver_net_events()
  // netcvode.cpp::NetCvode::check_thresh(NrnThread*)
  if (!benchmark) {
    if (branch->soma_) {
      // perform step-begin operations, eg send updates, wait for dependencies
      synchronizer_->StepSync(branch);

      // Soma waits for AIS to have threshold V value updated
      floble_t threshold_v;
      HinesSolver::SynchronizeThresholdV(branch, &threshold_v);
      if (branch->soma_->CheckAPthresholdAndTransmissionFlag(threshold_v))
        spikes_lco = branch->soma_->SendSpikes(nt->_t);
    } else if (branch->thvar_ptr_)
      // Axon Initial Segment send threshold  V to parent
      HinesSolver::SynchronizeThresholdV(branch);
  }

  // netcvode.cpp::NetCvode::deliver_net_events()
  t += .5 * dt;
  if (!benchmark) branch->DeliverEvents(t);  // events in 1st half-step
  branch->FixedPlayContinuous();
  SetupTreeMatrix(branch);
  HinesSolver::SolveTreeMatrix(branch);

  // update ions currents based on RHS and dI/dV
  second_order_cur(nt, input_params_->second_order_);

  ////// fadvance_core.c : update()
  HinesSolver::UpdateVoltagesWithRHS(branch);
  // TODO branch can be placed after the next operation

  // update capacitance currents based on RHS and dI/dV
  branch->CallModFunction(Mechanism::ModFunctions::kCurrentCapacitance);

  ////// fadvance_core.::nrn_fixed_step_lastpart()
  // callModFunction(Mechanism::ModFunction::jacob);
  t += .5 * dt;
  // TODO: commenting the call below changes nothing
  //(it changes the variables used in the current function only)
  branch->FixedPlayContinuous();
  branch->CallModFunction(Mechanism::ModFunctions::kState);
  branch->CallModFunction(Mechanism::ModFunctions::kAfterSolve);
  branch->CallModFunction(Mechanism::ModFunctions::kBeforeStep);
  if (!benchmark) branch->DeliverEvents(t);  // events in 2nd half-step

  return spikes_lco;
}

void BackwardEuler::RunFunction2(hpx_t lco, const Branch *b, const int func_id)
{
    //call thread locally
    hpx_call(HPX_HERE, BackwardEuler::RunFunction, lco,
             &b, sizeof(Branch*), &func_id, sizeof(func_id));
}


hpx_action_t BackwardEuler::RunFunction = 0;
int BackwardEuler::RunFunction_handler(const int nargs, const void *args[],
                                       const size_t[]) {
NEUROX_MEM_PIN(uint64_t);
assert(nargs==2);
Branch *branch = *(Branch**) args[0];
const int process_id = *(const int*) args[1];

switch (process_id)
{
  case 3:
    HinesSolver::SetupMatrixRHS(branch);
    break;
  case 4:
    ///branch->CallModFunction(Mechanism::ModFunctions::kJacob); //not used
    branch->CallModFunction(Mechanism::ModFunctions::kJacobCapacitance);
    HinesSolver::SetupMatrixDiagonal(branch);
    break;
  case 5:
    second_order_cur(branch->nt_, input_params_->second_order_);
    break;
  case 6:
    HinesSolver::UpdateVoltagesWithRHS(branch);
    break;
  case 7:
    branch->CallModFunction(Mechanism::ModFunctions::kCurrentCapacitance);
    break;
  default:
    assert(0);
  break;
}
NEUROX_MEM_UNPIN;
}

hpx_t BackwardEuler::StepParallel(Branch *branch) {
    NrnThread *nt = branch->nt_;
    double &t = nt->_t;
    const double dt = nt->_dt;
    hpx_t spikes_lco = HPX_NULL;
    hpx_t run_lco = HPX_NULL;

    // netcvode.cpp::NetCvode::deliver_net_events()
    //updates ml->data for point-processes
    branch->DeliverEvents(t + .5*dt);  // events in 1st half-step

    // sets vecplay->*pd which is ml->pdata for stimulus
    branch->FixedPlayContinuous( t + .5*dt);

    if (branch->soma_) {
        // perform step-begin operations, eg send updates, wait for dependencies
        synchronizer_->StepSync(branch);

        // Soma waits for AIS to have threshold V value updated
        floble_t threshold_v;
        HinesSolver::SynchronizeThresholdV(branch, &threshold_v);
        if (branch->soma_->CheckAPthresholdAndTransmissionFlag(threshold_v))
          spikes_lco = branch->soma_->SendSpikes(nt->_t);
    } else if (branch->thvar_ptr_)
        // Axon Initial Segment send threshold  V to parent
        HinesSolver::SynchronizeThresholdV(branch);

    t += .5 *dt;

    // zeroes RHS and D
    HinesSolver::ResetArray(branch, nt->_actual_rhs);
    HinesSolver::ResetArray(branch, nt->_actual_d);

    // zeroes and updates RHS and D from V and mech state
    //branch->CallModFunction(Mechanism::ModFunctions::kBeforeBreakpoint);
    branch->CallModFunction(Mechanism::ModFunctions::kCurrent);

    // 3: updates RHS from A B and V
    // HinesSolver::SetupMatrixRHS(branch);

    // 4: updates D from all capacitors mechanisms stats
    // branch->CallModFunction(Mechanism::ModFunctions::kJacob); //not used
    // branch->CallModFunction(Mechanism::ModFunctions::kJacobCapacitance);
    // 4.1 :updates D from A and B
    // HinesSolver::SetupMatrixDiagonal(branch);

    run_lco = hpx_lco_and_new(2);
    RunFunction2(run_lco, branch, 3);
    RunFunction2(run_lco, branch, 4);
    hpx_lco_wait_reset(run_lco);

    // calculates final RHS from D //requires all values to be available!
    HinesSolver::SolveTreeMatrix(branch);

    //5: update ions currents based on RHS and dI/dV
    //second_order_cur(branch->nt_, input_params_->second_order_);

    //6: updates V based on RHS
    //HinesSolver::UpdateVoltagesWithRHS(branch);

    //7: update capacitance currents based on RHS and dI/dV
    //branch->CallModFunction(Mechanism::ModFunctions::kCurrentCapacitance);

    run_lco = hpx_lco_and_new(3);
    RunFunction2(run_lco, branch, 5);
    RunFunction2(run_lco, branch, 6);
    RunFunction2(run_lco, branch, 7);
    hpx_lco_wait_reset(run_lco);

    t += .5 * dt;
    // TODO: commenting the call below changes nothing
    // (it changes the variables used in the current function only)
    // branch->FixedPlayContinuous();
    branch->CallModFunction(Mechanism::ModFunctions::kState);
    //branch->CallModFunction(Mechanism::ModFunctions::kAfterSolve);
    //branch->CallModFunction(Mechanism::ModFunctions::kBeforeStep);
    branch->DeliverEvents(t);  // events in 2nd half-step

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

void BackwardEuler::RegisterHpxActions() {
  wrappers::RegisterMultipleVarAction(RunFunction, RunFunction_handler);
}
