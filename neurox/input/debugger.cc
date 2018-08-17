#include "neurox/neurox.h"

#include <stdio.h>
#include <algorithm>
#include <list>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <tuple>

#include "coreneuron/nrnconf.h"
#include "coreneuron/nrniv/ivocvect.h"
#include "coreneuron/nrniv/netcon.h"
#include "coreneuron/nrniv/netcvode.h"
#include "coreneuron/nrniv/nrn_assert.h"
#include "coreneuron/nrniv/nrn_setup.h"
#include "coreneuron/nrniv/nrniv_decl.h"
#include "coreneuron/nrniv/vrecitem.h"
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrnoc/nrnoc_decl.h"
#include "coreneuron/nrnoc/nrnoc_decl.h"  //nrn_is_ion()
#include "coreneuron/utils/memory_utils.h"

using namespace std;
using namespace neurox::input;
using namespace neurox::synchronizers;
using namespace neurox::interpolators;

bool Debugger::IsEqual(floble_t a, floble_t b, bool roughlyEqual) {
  const double epsilon = input_params_->graph_mechs_parallelism_ ? 5e-4 : 1e-8;
  // NOTE: current functions for mechs SK_E2 (144) and Nap_Et2 (123) depends on
  // V that change opening vars that change V and so on... after few steps this
  // causes high difference in results;
  return roughlyEqual ? fabs(a - b) < epsilon : a == b;
}

void Debugger::CompareMechanismsFunctions() {
#if !defined(NDEBUG)
  DebugMessage(
      "neurox::Input::CoreNeuron::Debugger::CompareMechanismsFunctions...\n");
  for (int m = 0; m < neurox::mechanisms_count_; m++) {
    int index = mechanisms_[m]->type_;
    Memb_func &mf_cn = memb_func[index];            // coreneuron
    Memb_func &mf_nx = mechanisms_[m]->memb_func_;  // neurox

    // constructur, destructors are different
    // assert (mf_cn.alloc == mf_nx.alloc);
    // assert (mf_cn.destructor == mf_nx.destructor);

    // file reader is not used (we read directly from Coreneuron)
    // assert (mechanisms[m]->nrn_bbcore_read == nrn_bbcore_read_[index]);

    // cap calls MechFunctions::jacobCapacitance and currentCapacitance instead
    if (index != CAP) {
      assert(mf_cn.current == mf_nx.current);
      assert(mf_cn.jacob == mf_nx.jacob);
    }

    assert(mf_cn.state == mf_nx.state);
    assert(mf_cn.initialize == mf_nx.initialize);

    // If fails below, probably forgot to disable threading table when
    // generating circuit
    assert(mf_cn.setdata_ == mf_nx.setdata_);
    assert(mf_cn.thread_cleanup_ == mf_nx.thread_cleanup_);
    assert(mf_cn.thread_mem_init_ == mf_nx.thread_mem_init_);
    assert(mf_cn.thread_size_ == mf_nx.thread_size_);
    assert(mf_cn.thread_table_check_ == mf_nx.thread_table_check_);

    for (int i = 0; i < BEFORE_AFTER_SIZE; i++)
      if (nrn_threads[0].tbl[i])
        assert(mechanisms_[m]->before_after_functions_[i] ==
               nrn_threads[0].tbl[i]->bam->f);
  }
#endif
}

void Debugger::StepAfterStepFinitialize(Branch *b, NrnThread *nth) {
  b->CallModFunction(Mechanism::ModFunctions::kThreadTableCheck);
  b->InitVecPlayContinous();

  // coreneuron
  t = 0.;
  dt2thread(-1.);
  nrn_thread_table_check();
  clear_event_queue();
  nrn_spike_exchange_init();
  nrn_play_init();

  /****************/ CompareBranch2(b); /*****************/

  b->DeliverEvents(b->nt_->_t);
  for (int n = 0; n < b->nt_->end; n++)
    b->nt_->_actual_v[n] = input_params_->voltage_;

  // coreneuron
  nrn_deliver_events(nth); /* The play events at t=0 */

  for (int i = 0; i < nth->end; ++i)
    nth->_actual_v[i] = input_params_->voltage_;

  /****************/ CompareBranch2(b); /*****************/

  b->CallModFunction(Mechanism::ModFunctions::kBeforeInitialize);

  nrn_ba(nth, BEFORE_INITIAL);

  /****************/ CompareBranch2(b); /*****************/

  b->CallModFunction(Mechanism::ModFunctions::kInitialize);

  NrnThreadMembList *tml;
  for (tml = nth->tml; tml; tml = tml->next) {
    mod_f_t s = memb_func[tml->index].initialize;
    if (s) {
      (*s)(nth, tml->ml, tml->index);
    }
  }
  init_net_events();

  /****************/ CompareBranch2(b); /*****************/

  b->CallModFunction(Mechanism::ModFunctions::kAfterInitialize);
  b->DeliverEvents(b->nt_->_t);

  // coreneuron

  nrn_ba(nth, AFTER_INITIAL);
  nrn_deliver_events(nth); /* The INITIAL sent events at t=0 */

  /****************/ CompareBranch2(b); /*****************/

  BackwardEuler::SetupTreeMatrix(b);

  setup_tree_matrix_minimal(nth);

  /****************/ CompareBranch2(b); /*****************/

  b->CallModFunction(Mechanism::ModFunctions::kBeforeStep);
  b->DeliverEvents(b->nt_->_t);

  nrn_ba(nth, BEFORE_STEP);
  nrn_deliver_events(nth);

  /****************/ CompareBranch2(b); /*****************/
}

void Debugger::StepAfterStepBackwardEuler(Branch *b, NrnThread *nth,
                                          int secondorder) {
  double dt = b->nt_->_dt;
  if (b->soma_) {
    synchronizer_->StepSync(b);
    // Soma waits for AIS to have threshold V value updated
    floble_t thresholdV = HinesSolver::GetAxonInitialSegmentVoltage(b);
    if (b->soma_->CheckAPthresholdAndTransmissionFlag(thresholdV))
      b->soma_->SendSpikes(b->nt_->_t);
    // TODO sendSpikes LCO must be waited
  } else if (b->thvar_ptr_)
    // Axon Initial Segment send threshold  V to parent
    HinesSolver::GetAxonInitialSegmentVoltage(b);

  b->nt_->_t += .5 * dt;
  b->DeliverEvents(b->nt_->_t);

  // coreneuron
  assert(b->nt_->cj == nth->cj);
  dt2thread(dt);
  ;
  deliver_net_events(nth);
  nth->_t += .5 * nth->_dt;

  /****************/ CompareBranch2(b); /*****************/

  b->FixedPlayContinuous();
  BackwardEuler::SetupTreeMatrix(b);
  HinesSolver::SolveTreeMatrix(b);

  // coreneuron
  fixed_play_continuous(nth);
  setup_tree_matrix_minimal(nth);
  nrn_solve_minimal(nth);

  /****************/ CompareBranch2(b); /*****************/

  second_order_cur(b->nt_, input_params_->second_order_);
  floble_t secondOrderMultiplier = input_params_->second_order_ ? 2.0 : 1.0;
  for (int i = 0; i < b->nt_->end; i++)
    b->nt_->_actual_v[i] += secondOrderMultiplier * b->nt_->_actual_rhs[i];
  b->CallModFunction(Mechanism::ModFunctions::kCurrentCapacitance);

  second_order_cur(nth, secondorder);
  update(nth);

  /****************/ CompareBranch2(b); /*****************/

  b->nt_->_t += .5 * dt;
  b->FixedPlayContinuous();

  nth->_t += .5 * nth->_dt;  // nrn_fixed_step_lastpart(nth);
  fixed_play_continuous(nth);

  /****************/ CompareBranch2(b); /*****************/

  b->CallModFunction(Mechanism::ModFunctions::kState);

  nonvint(nth);

  /****************/ CompareBranch2(b); /*****************/

  b->CallModFunction(Mechanism::ModFunctions::kAfterSolve);
  b->CallModFunction(Mechanism::ModFunctions::kBeforeStep);
  b->DeliverEvents(b->nt_->_t);

  nrn_ba(nth, AFTER_SOLVE);
  nrn_ba(nth, BEFORE_STEP);
  nrn_deliver_events(nth);

  /****************/ CompareBranch2(b); /*****************/
}

// fadvance_core.c::fixed_step_minimal(...)
void Debugger::FixedStepMinimal2(NrnThread *nth, int secondorder) {
  // fadvance_core.c::dt2thread
  if (secondorder)
    nth->cj = 2.0 / nth->_dt;
  else
    nth->cj = 1.0 / nth->_dt;

  // fadvance_core.c::nrn_fixed_step_thread(...)
  deliver_net_events(nth);
  nth->_t += .5 * nth->_dt;
  if (nth->ncell) {
    fixed_play_continuous(nth);
    setup_tree_matrix_minimal(nth);
    nrn_solve_minimal(nth);
    second_order_cur(nth, secondorder);
    update(nth);
  }
  if (!nrn_have_gaps) {
    nrn_fixed_step_lastpart(nth);
  }
}

void Debugger::CompareAllBranches() {
#if !defined(NDEBUG)
  if (/*input_params_->branch_parallelism_ ||  */
      /* (branch parallelism is possible to be compared if there are no branches
       * and this usecase is handled by CoreNeuron::Debugger::CompareBranch) */
      input_params_->load_balancing_ ||
      input_params_->interpolator_ != InterpolatorIds::kBackwardEuler)
    return;

  DebugMessage("neurox::Input::CoreNeuron::Debugger::CompareBranch...\n");
  neurox::wrappers::CallAllNeurons(input::Debugger::CompareBranch);
#endif
}

void Debugger::CompareBranch2(Branch *branch) {
  assert(branch->soma_);  // only non-branched neurons

  int nrnThreadId = branch->nt_->id;
  assert(sizeof(floble_t) == sizeof(double));  // only works with doubles!
  NrnThread &nt = nrn_threads[nrnThreadId];

  bool multiMex = branch->mechs_graph_ != NULL;

  assert(IsEqual(branch->nt_->_t, nt._t, true));
  assert(branch->nt_->_ndata == nt._ndata);
  assert(secondorder == input_params_->second_order_);
  assert(branch->soma_->threshold_ == nt.presyns[0].threshold_);
  assert(branch->soma_->gid_ == nt.presyns[0].gid_);
  assert(IsEqual(*(branch->thvar_ptr_),
                 nt._actual_v[nt.presyns[0].thvar_index_], multiMex));

  // vecplay
  assert(branch->nt_->n_vecplay == nt.n_vecplay);
  for (int i = 0; i < nt.n_vecplay; i++) {
    void *vpx_ptr = branch->nt_->_vecplay[i];
    VecplayContinuousX *vpx = reinterpret_cast<VecplayContinuousX *>(vpx_ptr);

    PlayRecord *prc = (PlayRecord *)nt._vecplay[i];
    VecPlayContinuous *vpc = (VecPlayContinuous *)prc;

    assert(prc->ith_ == branch->nt_->id);  // single neuron per NrnThread
    assert(vpx->last_index_ == vpc->last_index_);
    assert(vpx->discon_index_ == vpc->discon_index_);
    assert(vpx->size_ == vpc->y_->size());
    assert(vpx->size_ == vpc->t_->size());
    assert(IsEqual(*(vpx->pd_), *(vpc->pd_), multiMex));
    assert((vpx->discon_indices_ == NULL) == (vpc->discon_indices_ == NULL));
    for (size_t j = 0; j < vpx->size_; j++) {
      assert(vpx->y_[j] == vpc->y_->data()[j]);  // constants
      assert(vpx->t_[j] == vpc->t_->data()[j]);  // constants
      if (vpx->discon_indices_ != NULL) {
        assert(vpx->discon_indices_[j] == vpc->discon_indices_->data()[j]);
      }  // constants
    }
  }

  // make sure morphology is correct
  assert(branch->nt_->end == nt.end);
  for (int i = 0; i < 6; i++)                   // RHS, D, A, B, V and area
    for (int j = 0; j < branch->nt_->end; j++)  // for all compartments
    {
      int offset = tools::Vectorizer::SizeOf(nt.end) * i + j;
      assert(IsEqual(nt._data[offset], branch->nt_->_data[offset], multiMex));
    }

  // dparam_semantics
  for (int m = 0; m < neurox::mechanisms_count_; m++) {
    int type = mechanisms_[m]->type_;
    for (int i = 0; i < mechanisms_[m]->pdata_size_; i++) {
      assert(mechanisms_[m]->memb_func_.dparam_semantics[i] ==
             memb_func[type].dparam_semantics[i]);
    }
  }

  for (offset_t i = 0; i < branch->nt_->end; i++) {
    assert(nt._actual_a[i] == branch->nt_->_actual_a[i]);        // constants
    assert(nt._actual_b[i] == branch->nt_->_actual_b[i]);        // constants
    assert(nt._actual_area[i] == branch->nt_->_actual_area[i]);  // constants
    assert(IsEqual(nt._actual_d[i], branch->nt_->_actual_d[i], multiMex));
    assert(IsEqual(nt._actual_v[i], branch->nt_->_actual_v[i], multiMex));
    assert(IsEqual(nt._actual_rhs[i], branch->nt_->_actual_rhs[i], multiMex));
    if (branch->nt_->_v_parent_index) {
      assert(nt._v_parent_index[i] == branch->nt_->_v_parent_index[i]);
    }
  }

  // make sure weights are correct
  assert(nt.n_weight == branch->nt_->n_weight);
  // for (int n=0; n< nt.n_weight; n++)
  //{   assert(nt.weights[n] == branch->nt->weights[n]); } //order is changed!

  // make sure netcons are correct
  size_t netconsCount = 0;
  for (auto nc_it : branch->netcons_) {
    int srcGid = nc_it.first;
    netconsCount += nc_it.second.size();
  }
  assert(nt.n_netcon == netconsCount);

  int vdataOffset = 0;
  for (NrnThreadMembList *tml = nt.tml; tml != NULL;
       tml = tml->next)  // For every mechanism
  {
    int type = tml->index;
    int m = mechanisms_map_[type];
    Memb_list *ml = tml->ml;  // Mechanisms application to each compartment
    Memb_list &instances = branch->mechs_instances_[m];
    assert(ml->nodecount == instances.nodecount);
    short dataSize = mechanisms_[m]->data_size_;
    short pdataSize = mechanisms_[m]->pdata_size_;
    assert(pdataSize == nrn_prop_dparam_size_[type]);
    assert(dataSize == nrn_prop_param_size_[type]);
    if (DataLoader::HardCodedMechanismHasNoInstances(tml->index)) {
      assert(ml->nodeindices == nullptr and instances.nodeindices == nullptr);
      assert(ml->nodecount == instances.nodecount);
    } else
      for (int n = 0; n < ml->nodecount; n++)  // for every mech instance
      {
        assert(ml->nodeindices[n] == instances.nodeindices[n]);
        for (int i = 0; i < dataSize; i++) {
#if LAYOUT == 1
          int offset = neurox::mechanisms_[m]->data_size_ * n + i;
#else
          int offset = tools::Vectorizer::SizeOf(ml->nodecount) * i + n;
#endif
          assert(IsEqual(ml->data[offset], instances.data[offset], multiMex));
        }

        for (int i = 0; i < pdataSize; i++) {
#if LAYOUT == 1
          int offset = pdataSize * n + i;
#else
          int offset = tools::Vectorizer::SizeOf(ml->nodecount) * i + n;
#endif
          int ptype = memb_func[type].dparam_semantics[i];
          bool isPointer = ptype == -1 || (ptype > 0 && ptype < 1000);
          if (isPointer) {
            assert(IsEqual(nt._data[ml->pdata[offset]],
                           branch->nt_->_data[instances.pdata[offset]],
                           multiMex));
          }
          assert(ml->pdata[offset] == instances.pdata[offset]);
        }

        /* We comment this because it runs for NULL presyn
        if (mechanisms[m]->pntMap)
            //compare point_processes (index 1)
            Point_process * pp = (Point_process *) nt._vdata[vdataOffset+1];
            Point_process * pp2 = (Point_process *)
        branch->vdata[vdataOffset+1]; assert(pp->_type == pp2->_type );
            assert(pp->_i_instance == pp2->_i_instance );
            assert(pp->_tid == pp2->_tid );
            assert(pp->_presyn == pp2->_presyn );
            vdataOffset+= mechanisms[m]->vdataSize;
        }
        */
      }
  }
}

hpx_action_t Debugger::FixedStepMinimal = 0;
int Debugger::FixedStepMinimal_handler(const int *steps_ptr, const size_t) {
  NEUROX_MEM_PIN(uint64_t);
  for (int n = 0; n < nrn_nthread; n++)
    for (int i = 0; i < *steps_ptr; i++)
      Debugger::FixedStepMinimal2(&nrn_threads[n],
                                  input_params_->second_order_);
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t Debugger::Finitialize = 0;
int Debugger::Finitialize_handler() {
  NEUROX_MEM_PIN(uint64_t);
  if (input_params_->interpolator_ == InterpolatorIds::kBackwardEuler)
    nrn_finitialize(input_params_->voltage_ != 1000., input_params_->voltage_);
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t Debugger::ThreadTableCheck = 0;
int Debugger::ThreadTableCheck_handler() {
  // beginning of fadvance_core.c::nrn_fixed_step_group_minimal
  NEUROX_MEM_PIN(uint64_t);
  if (input_params_->interpolator_ == InterpolatorIds::kBackwardEuler) {
    dt2thread(dt);  // does nothing
    nrn_thread_table_check();
  }
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t Debugger::NrnSpikeExchange = 0;
int Debugger::NrnSpikeExchange_handler() {
  NEUROX_MEM_PIN(uint64_t);
#if NRNMPI == 1
  if (input_params_->interpolator_ == InterpolatorIds::kBackwardEuler)
    nrn_spike_exchange(nrn_threads);
#endif
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t Debugger::CompareBranch = 0;
int Debugger::CompareBranch_handler() {
  NEUROX_MEM_PIN(Branch);
  if ((input_params_->branch_parallelism_ && local->branch_tree_ &&
       local->branch_tree_->branches_count_ > 0) ||
      input_params_->load_balancing_ ||
      input_params_->interpolator_ != InterpolatorIds::kBackwardEuler)
    return neurox::wrappers::MemoryUnpin(target);
  CompareBranch2(local);  // not implemented for branch-parallelism
  return neurox::wrappers::MemoryUnpin(target);
}

void Debugger::RunCoreneuronAndCompareAllBranches() {
#if !defined(NDEBUG)
  if (input_params_->branch_parallelism_ || input_params_->load_balancing_ ||
      input_params_->interpolator_ != InterpolatorIds::kBackwardEuler)
    return;

  // parallel execution only (serial execs are compared on-the-fly
  if (neurox::ParallelExecution()) {
    const int total_steps = BackwardEuler::GetTotalSteps();
    const int comm_steps = BackwardEuler::GetMinSynapticDelaySteps();
    DebugMessage(
        "neurox::re-running simulation in Coreneuron to compare final "
        "result...\n");
    for (int s = 0; s < total_steps; s += comm_steps) {
      hpx_bcast_rsync(neurox::input::Debugger::FixedStepMinimal, &comm_steps,
                      sizeof(int));
      hpx_bcast_rsync(neurox::input::Debugger::NrnSpikeExchange);
    }
    neurox::input::Debugger::CompareAllBranches();
  }
#endif
}

void Debugger::SingleNeuronStepAndCompare(NrnThread *nt, Branch *b,
                                          char secondorder) {
#if !defined(NDEBUG)
  if (input_params_->branch_parallelism_ || input_params_->load_balancing_ ||
      neurox::ParallelExecution())
    return;  // can't be compared

  input::Debugger::FixedStepMinimal2(nt, secondorder);
  input::Debugger::CompareBranch2(b);
#endif
}

void Debugger::RegisterHpxActions() {
  wrappers::RegisterZeroVarAction(Debugger::CompareBranch,
                                  Debugger::CompareBranch_handler);
  wrappers::RegisterZeroVarAction(Debugger::Finitialize,
                                  Debugger::Finitialize_handler);
  wrappers::RegisterZeroVarAction(Debugger::NrnSpikeExchange,
                                  Debugger::NrnSpikeExchange_handler);
  wrappers::RegisterZeroVarAction(Debugger::ThreadTableCheck,
                                  Debugger::ThreadTableCheck_handler);
  wrappers::RegisterSingleVarAction<int>(Debugger::FixedStepMinimal,
                                         Debugger::FixedStepMinimal_handler);
  //  NEUROX_REGISTER_ACTION(NEUROX_ACTION_ZERO_VAR, Debugger::CompareBranch);
  //  NEUROX_REGISTER_ACTION(NEUROX_ACTION_SINGLE_VAR,
  //  Debugger::FixedStepMinimal);
}
