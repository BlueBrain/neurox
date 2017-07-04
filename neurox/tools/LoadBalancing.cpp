#include "neurox/neurox.h"

using namespace std;
using namespace neurox;

double * Tools::LoadBalancing::loadBalancingTable = nullptr;
hpx_t Tools::LoadBalancing::loadBalancingMutex = HPX_NULL;

hpx_action_t Tools::LoadBalancing::queryLoadBalancingTable = 0;
int Tools::LoadBalancing::queryLoadBalancingTable_handler(const int nargs, const void *args[], const size_t[])
{
    /**
     * nargs=1 or 2 where
     * args[0] = elapsedTime
     * args[1] = rank (if any)
     */
    neurox_hpx_pin(uint64_t);
    assert(nargs==1 || nargs==2);
    assert(hpx_get_my_rank()==0); //only one loadBalancingTable and only in rank zero
    const double elapsedTime = *(const double*)args[0];

    if (nargs==2) //this neuron already has a rank allocated, update it's entry
    {
        const int rank = *(const int*)args[1];
        hpx_lco_sema_p(loadBalancingMutex);
        loadBalancingTable[rank] += elapsedTime;
        hpx_lco_sema_v_sync(loadBalancingMutex);
        neurox_hpx_unpin;
    }
    else
    {
        double minElapsedTime=99999999999;
        int rank=-1;
        hpx_lco_sema_p(loadBalancingMutex);
        for (int r=0; r<hpx_get_num_ranks(); r++)
            if (loadBalancingTable[r]<minElapsedTime)
            {
                minElapsedTime = loadBalancingTable[r];
                rank=r;
            }
        loadBalancingTable[rank] += elapsedTime;
        hpx_lco_sema_v_sync(loadBalancingMutex);
        neurox_hpx_unpin_continue(rank);
    }
    neurox_hpx_unpin;
}

Tools::LoadBalancing::LoadBalancing()
{
    if (hpx_get_my_rank()==0)
    {
        loadBalancingMutex = hpx_lco_sema_new(1);
        loadBalancingTable = new double[hpx_get_num_ranks()];
        for (int r=0; r<hpx_get_num_ranks(); r++)
            loadBalancingTable[r]=0;
    }
}

Tools::LoadBalancing::~LoadBalancing()
{
    hpx_lco_delete_sync(loadBalancingMutex);
    delete [] loadBalancingTable;
}

void Tools::LoadBalancing::print()
{
    for (int r=0; r<hpx_get_num_ranks(); r++)
        printf("-- rank %d : %.6f ms\n", r, loadBalancingTable[r]);
}

void Tools::LoadBalancing::registerHpxActions()
{
    neurox_hpx_register_action(neurox_several_vars_action, queryLoadBalancingTable);
}
