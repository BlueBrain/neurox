#include "neurox/neurox.h"

using namespace std;
using namespace neurox;

double *tools::LoadBalancing::load_balancing_table_ = nullptr;
hpx_t tools::LoadBalancing::load_balancing_mutex_ = HPX_NULL;
double tools::LoadBalancing::total_mech_instances_runtime_ = 0;

hpx_action_t tools::LoadBalancing::QueryLoadBalancingTable = 0;
int tools::LoadBalancing::QueryLoadBalancingTable_handler(const int nargs,
                                                          const void *args[],
                                                          const size_t[]) {
  /**
   * nargs=1 or 2 where
   * args[0] = elapsedTime
   * args[1] = rank (if any)
   */
  NEUROX_MEM_PIN(uint64_t);
  assert(nargs == 1 || nargs == 2);

  // only one loadBalancingTable and only in rank zero
  assert(hpx_get_my_rank() == 0);
  const double elapsed_time = *(const double *)args[0];

  // this neuron already has a rank allocated, update it's entry
  if (nargs == 2) {
    const int rank = *(const int *)args[1];
    hpx_lco_sema_p(load_balancing_mutex_);
    load_balancing_table_[rank] += elapsed_time;
    hpx_lco_sema_v_sync(load_balancing_mutex_);
    return wrappers::MemoryUnpin(target);
  } else {
    double min_elapsed_time = 99999999999;
    int rank = -1;
    hpx_lco_sema_p(load_balancing_mutex_);
    for (int r = 0; r < hpx_get_num_ranks(); r++)
      if (load_balancing_table_[r] < min_elapsed_time) {
        min_elapsed_time = load_balancing_table_[r];
        rank = r;
      }
    load_balancing_table_[rank] += elapsed_time;
    hpx_lco_sema_v_sync(load_balancing_mutex_);
    return wrappers::MemoryUnpin(target, rank);
  }
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

  printf("neurox::tools::LoadBalancing::PrintTable()\n");
  for (int r = 0; r < hpx_get_num_ranks(); r++)
    printf("- rank %d : %.6f ms\n", r, load_balancing_table_[r]);
}

double tools::LoadBalancing::GetWorkPerBranchSubsection(
    const double neuron_time, const int neurons_count) {
  // for a single-thread CPU with a single neuron..
  double work_per_section = neuron_time / neurons_count;

  // for a multi-thread CPU with a single neuron...
  work_per_section /= wrappers::NumThreads();

  // for a multi-thread CPU with several neurons...
  work_per_section *= neurons_count;

  // if graph parallelism is available, increase the work by a factor of...?
  if (input_params_->graph_mechs_parallelism_) work_per_section *= 10;
}

void tools::LoadBalancing::AddToTotalMechInstancesRuntime(double runtime)
{
    total_mech_instances_runtime_ += runtime;
}

double tools::LoadBalancing::GetWorkloadPerMechInstancesThread()
{
    return total_mech_instances_runtime_ * kMechInstancesPercentagePerComputeUnit;
}

void tools::LoadBalancing::RegisterHpxActions() {
  wrappers::RegisterMultipleVarAction(QueryLoadBalancingTable,
                                      QueryLoadBalancingTable_handler);
}
