#include "neurox/neurox.h"
#include <cstring>
#include <algorithm>

Branch::Branch(const int n, const double *a, const double *b, const double *d,
               const double *v, const double *rhs, const double *area,
               const int m, const int * mechsCounts, const double *data,
               const Datum *pdata, const int childrenCount, const hpx_t * children)
    :n(n), m(m), childrenCount(childrenCount)
{
    this->a = new double[n];
    this->b = new double[n];
    this->d = new double[n];
    this->v = new double[n];
    this->rhs = new double[n];
    this->area = new double[n];
    this->mechsOffsets = new int[brain->mechsTypesCount];
    this->children = new hpx_t[childrenCount];

    memcpy(this->a,a,n*sizeof(double));
    memcpy(this->b,b,n*sizeof(double));
    memcpy(this->d,d,n*sizeof(double));
    memcpy(this->v,v,n*sizeof(double));
    memcpy(this->rhs,rhs,n*sizeof(double));
    memcpy(this->area,area,n*sizeof(double));
    memcpy(this->children, children, childrenCount*sizeof(hpx_t));

    //calculate offsets based on count
    int dataSize=0, pdataSize=0;
    for (int i=0; i<brain->mechsTypesCount; i++)
    {
        dataSize += brain->mechsTypes[i].dataSize * mechsCounts[i];
        pdataSize += brain->mechsTypes[i].pdataSize * mechsCounts[i];
        mechsOffsets[i] = i==0 ? 0 : mechsOffsets[i-1]+ mechsCounts[i];
    }

    this->data = new double[dataSize];
    this->pdata = new Datum[pdataSize];
    memcpy(this->data, data, dataSize*sizeof(double));
    memcpy(this->pdata, pdata, pdataSize*sizeof(Datum));
}

Branch::~Branch()
{
    delete [] b;
    delete [] d;
    delete [] a;
    delete [] v;
    delete [] rhs;
    delete [] area;
    delete [] mechsOffsets;
    delete [] data;
    delete [] pdata;
    delete [] children;
}

hpx_action_t Branch::init = 0;
int Branch::init_handler(const int n, const double *a, const double *b, const double *d,
                               const double *v, const double *rhs, const double *area,
                               const int m, const int * mechsCount, const double *data,
                               const Datum *pdata, const int childrenCount, const hpx_t * children)
{
    //Make sure message arrived correctly, and pin memory
    hpx_t branch_addr = hpx_thread_current_target();
    Branch * branch = NULL;
    if (!hpx_gas_try_pin(branch_addr, (void**) &branch))
        return HPX_RESEND;

    //do the work
    branch = new Branch(n,a,b,d,v,rhs,area,m,mechsCount, data, pdata, childrenCount, children);

    //unpin and return success
    hpx_gas_unpin(branch_addr);
    return HPX_SUCCESS;
}

hpx_action_t Branch::finitialize = 0;
int Branch::finitialize_handler()
{
    //Make sure message arrived correctly, and pin memory
    hpx_t branch_addr = hpx_thread_current_target();
    Branch * branch = NULL;
    if (!hpx_gas_try_pin(branch_addr, (void**) &branch))
        return HPX_RESEND;

    //set by finitialize.c:nrn_finitialize()
    for (int i=0; i<branch->n; i++)
        branch->v[i]=inputParams->voltage;

    setupTreeMatrixMinimal(branch);

    //call it for every sub-branch
    hpx_par_for_sync(
       [&] (int i, void*) { hpx_call_sync(branch->children[i], Branch::finitialize);},
       0, branch->childrenCount, NULL);

    //unpin and return success
    hpx_gas_unpin(branch_addr);
    return HPX_SUCCESS;
}

void Branch::setupTreeMatrixMinimal(Branch * branch)
{
    //RIGHT HAND SIDE: treesetcore.c::nrn_rhs
    for (int i=0; i<branch->n; i++)
    {
        branch->rhs[i]=0;
        branch->d[i]=0;
    }

    //call current method in all mechanisms (inside the mod file it iterates across all applications)
    for (int i=0; i<brain->mechsTypesCount; i++)
        if (brain->mechsTypes[i].current)
            (brain->mechsTypes[i].current)(/*NrnThread nt, Memb_list ml. i*/);
             //TODO pass pdata/data offsets and count of mechs applications instead

    //TODO run nrn_ba(_nt, BEFORE_BREAKPOINT);

    //the internal axial currents: rhs += ai_j*(vi_j - vi)
    //The extracellular mechanism contribution is already done.
    for (int i=0; i<branch->n; i++)
    {
       //TODO
    }

    //TODO ---------------------------

    /* LEFT HAND SIDE: treesetcore.c::nrn_lhs
         cm*dvm/dt = -i(vm) + is(vi) + ai_j*(vi_j - vi)
         cx*dvx/dt - cm*dvm/dt = -gx*(vx - ex) + i(vm) + ax_j*(vx_j - vx)
       with a matrix so that the solution is of the form [dvm+dvx,dvx] on the right
       hand side after solving.
       This is a common operation for fixed step, cvode, and daspk methods
    */

    //call jacob method in all mechanisms
    for (int i=0; i<brain->mechsTypesCount; i++)
        if (brain->mechsTypes[i].jacob)
            (brain->mechsTypes[i].jacob)(/*NrnThread nt, Memb_list ml. i*/);

    //now the cap current can be computed because any change to cm by another model has taken effect
    //TODO nrn_cap_jaco

    //now add the axial currents
    for (int i=0; i<branch->n; i++)
    {

    }
}

void Branch::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  init, init_handler, HPX_INT, HPX_POINTER,
                        HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_INT, HPX_POINTER,
                        HPX_POINTER, HPX_POINTER, HPX_INT, HPX_POINTER);

    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  finitialize, finitialize_handler);
}
