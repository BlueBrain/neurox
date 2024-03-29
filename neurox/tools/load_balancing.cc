/*
# =============================================================================
# Copyright (c) 2015 - 2021-2022 Blue Brain Project/EPFL
#
# See top-level LICENSE file for details.
# =============================================================================.
*/
#include "neurox/neurox.h"

using namespace std;
using namespace neurox;

double *tools::LoadBalancing::load_balancing_table_ = nullptr;
hpx_t tools::LoadBalancing::load_balancing_mutex_ = HPX_NULL;

hpx_action_t tools::LoadBalancing::QueryLoadBalancingTable = 0;
int tools::LoadBalancing::QueryLoadBalancingTable_handler() {
  NEUROX_MEM_PIN(uint64_t);

  // only one loadBalancingTable and only in rank zero
  assert(hpx_get_my_rank() == 0);
  double min_elapsed_time = 99999999999;
  int rank = -1;
  hpx_lco_sema_p(load_balancing_mutex_);

  // return rank of least busy locality
  for (int r = 0; r < hpx_get_num_ranks(); r++)
    if (load_balancing_table_[r] < min_elapsed_time) {
      min_elapsed_time = load_balancing_table_[r];
      rank = r;
    }
  hpx_lco_sema_v_sync(load_balancing_mutex_);
  NEUROX_MEM_UNPIN_CONTINUE(rank);
}

hpx_action_t tools::LoadBalancing::UpdateLoadBalancingTable = 0;
int tools::LoadBalancing::UpdateLoadBalancingTable_handler(const int nargs,
                                                           const void *args[],
                                                           const size_t[]) {
  /**
   * nargs=2 where
   * args[0] = elapsed_time
   * args[1] = rank
   */
  NEUROX_MEM_PIN(uint64_t);
  assert(nargs == 2);

  // only one loadBalancingTable and only in rank zero
  assert(hpx_get_my_rank() == 0);
  const double elapsed_time = *(const double *)args[0];

  // this neuron already has a rank allocated, update it's entry
  const int rank = *(const int *)args[1];
  hpx_lco_sema_p(load_balancing_mutex_);
  load_balancing_table_[rank] += elapsed_time;
  hpx_lco_sema_v_sync(load_balancing_mutex_);
  return wrappers::MemoryUnpin(target);
  NEUROX_MEM_UNPIN;
}

tools::LoadBalancing::LoadBalancing() {
  if (hpx_get_my_rank() == 0) {
    load_balancing_mutex_ = hpx_lco_sema_new(1);
    load_balancing_table_ = new double[hpx_get_num_ranks()];
    for (int r = 0; r < hpx_get_num_ranks(); r++) load_balancing_table_[r] = 0;
  }
}

tools::LoadBalancing::~LoadBalancing() {
  hpx_lco_delete_sync(load_balancing_mutex_);
  delete[] load_balancing_table_;
}

void tools::LoadBalancing::PrintLoadBalancingTable() {
  if (load_balancing_table_ == nullptr) return;
  if (hpx_get_my_rank() != 0) return;

  printf("neurox::tools::LoadBalancing::PrintLoadBalancingTable:\n");
  for (int r = 0; r < hpx_get_num_ranks(); r++)
    printf("- rank %d : %.6f ms\n", r, load_balancing_table_[r]);
}

// TODO should be avg_neuron_time, not neuron_time (same in next func)!
double tools::LoadBalancing::GetMaxWorkPerBranchSubTree(
    const double neuron_time, const int my_neurons_count) {
  // intentionally made simillar to NEURON (see multisplit.hoc)
  return neuron_time / wrappers::NumThreads() *
         input_params_->subtree_complexity;
}

double tools::LoadBalancing::GetMaxWorkPerBranchSubSection(
    const double neuron_time, const int my_neurons_count) {
  if (input_params_->subsection_complexity == 0)
    return 0;  // 0 means "keep it as single subsection"

  // get an estimation of total workload in this locality
  double total_work = neuron_time * (double)my_neurons_count;

  // average the total work in this locality by all localityes
  double avg_locality_work = total_work / (double)wrappers::NumRanks();

  // scale to constant 1/k' (see paper)
  return avg_locality_work / input_params_->subsection_complexity;
}

double tools::LoadBalancing::GetWorkloadPerMechInstancesBlock(
    const double total_mech_instances_runtime) {
  return total_mech_instances_runtime *
         input_params_->mech_instance_percent_per_block;
}

void tools::LoadBalancing::RegisterHpxActions() {
  wrappers::RegisterZeroVarAction(QueryLoadBalancingTable,
                                  QueryLoadBalancingTable_handler);
  wrappers::RegisterMultipleVarAction(UpdateLoadBalancingTable,
                                      UpdateLoadBalancingTable_handler);
}
