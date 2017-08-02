#include "neurox/neurox.h"

using namespace std;
using namespace neurox;

double *tools::LoadBalancing::loadBalancingTable = nullptr;
hpx_t tools::LoadBalancing::loadBalancingMutex = HPX_NULL;

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
  assert(hpx_get_my_rank() ==
         0);  // only one loadBalancingTable and only in rank zero
  const double elapsedTime = *(const double *)args[0];

  if (nargs ==
      2)  // this neuron already has a rank allocated, update it's entry
  {
    const int rank = *(const int *)args[1];
    hpx_lco_sema_p(loadBalancingMutex);
    loadBalancingTable[rank] += elapsedTime;
    hpx_lco_sema_v_sync(loadBalancingMutex);
    NEUROX_MEM_UNPIN;
  } else {
    double minElapsedTime = 99999999999;
    int rank = -1;
    hpx_lco_sema_p(loadBalancingMutex);
    for (int r = 0; r < hpx_get_num_ranks(); r++)
      if (loadBalancingTable[r] < minElapsedTime) {
        minElapsedTime = loadBalancingTable[r];
        rank = r;
      }
    loadBalancingTable[rank] += elapsedTime;
    hpx_lco_sema_v_sync(loadBalancingMutex);
    NEUROX_MEM_UNPIN_CONTINUE(rank);
  }
  NEUROX_MEM_UNPIN;
}

tools::LoadBalancing::LoadBalancing() {
  if (hpx_get_my_rank() == 0) {
    loadBalancingMutex = hpx_lco_sema_new(1);
    loadBalancingTable = new double[hpx_get_num_ranks()];
    for (int r = 0; r < hpx_get_num_ranks(); r++) loadBalancingTable[r] = 0;
  }
}

tools::LoadBalancing::~LoadBalancing() {
  hpx_lco_delete_sync(loadBalancingMutex);
  delete[] loadBalancingTable;
}

void tools::LoadBalancing::PrintTable() {
  if (loadBalancingTable == nullptr) return;
  if (hpx_get_my_rank() != 0) return;

  printf("neurox::tools::LoadBalancing::PrintTable()\n");
  for (int r = 0; r < hpx_get_num_ranks(); r++)
    printf("- rank %d : %.6f ms\n", r, loadBalancingTable[r]);
}

void tools::LoadBalancing::RegisterHpxActions() {
  wrappers::RegisterMultipleVarAction(QueryLoadBalancingTable,
                                      QueryLoadBalancingTable_handler);
}
