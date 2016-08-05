#include "neurox/Neurox.h"
#include <cstring>
#include <algorithm>
#include <numeric>

using namespace NeuroX;
using namespace NeuroX::Solver;

hpx_action_t HinesSolver::setupMatrixRHS = 0;
int HinesSolver::setupMatrixRHS_handler(const double * parentV_ptr, const size_t)
{
    neurox_hpx_pin(Branch);
    const int n = local->n;
    double *a   = local->a;
    double *b   = local->b;
    double *v   = local->v;
    double *rhs = local->rhs;
    int branchesCount = local->branchesCount;

    double returnValue = -1; //contribution to upper branch
    double dv=-1;

    //for (i = i2; i < i3; ++i))
    if (!local->isSoma)
    {
        dv = *parentV_ptr-v[0];
        rhs[0] -= b[0]*dv;
        returnValue = a[0]*dv;
    }
    for (int i=1; i<local->n; i++)
    {
        dv = v[i-1]-v[i];
        rhs[i] -= b[i]*dv;
        rhs[i-1] += a[i]*dv;
    }

    if (branchesCount > 0)
    {
      //send/receive contribution to/from branches
      hpx_t * futures = new hpx_t[branchesCount];
      void  ** addrs  = new void*[branchesCount];
      size_t * sizes  = new size_t[branchesCount];
      double * values = new double[branchesCount];
      for (int c = 0; c < local->branchesCount; c++)
      {
        futures[c] = hpx_lco_future_new(sizeof(double));
        addrs[c]   = &values[c];
        sizes[c]   = sizeof(double);
        hpx_call(local->branches[c], HinesSolver::setupMatrixRHS, futures[c],
                 &v[n-1], sizeof(v[n-1]));
      }
      hpx_lco_get_all(branchesCount, futures, sizes, addrs, NULL);

      //received contributions, can now update value
      for (int c = 0; c < branchesCount; c++)
      {
        rhs[n-1] += values[c];
        hpx_lco_delete(futures[c], HPX_NULL);
      }

      delete [] futures;
      delete [] addrs;
      delete [] sizes;
      delete [] values;
    }

    if (!local->isSoma)
        neurox_hpx_unpin_continue(returnValue);
    neurox_hpx_unpin;
}

struct BackTriangFuture
{
    double rhs;
    double b;
}; ///> future value of the back-triangulation method

hpx_action_t HinesSolver::gaussianBackTriangulation = 0;
int HinesSolver::gaussianBackTriangulation_handler()
{
    neurox_hpx_pin(Branch);
    int n = local->n;
    double *a   = local->a;
    double *b   = local->b;
    double *d   = local->d;
    double *rhs = local->rhs;
    int branchesCount = local->branchesCount;

    hpx_t * futures = branchesCount ? new hpx_t[branchesCount]  : nullptr;
    void  ** addrs  = branchesCount ? new void*[branchesCount]  : nullptr;
    size_t * sizes  = branchesCount ? new size_t[branchesCount] : nullptr;

    BackTriangFuture * values = branchesCount ? new BackTriangFuture[branchesCount] : nullptr;

    for (int c = 0; c < branchesCount; c++)
    {
        futures[c] = hpx_lco_future_new(sizeof (BackTriangFuture));
        addrs[c]   = &values[c];
        sizes[c]   = sizeof(BackTriangFuture);
        hpx_call(local->branches[c], HinesSolver::gaussianBackTriangulation, futures[c]);
    }

    if (branchesCount > 0) //required or fails
        hpx_lco_get_all(branchesCount, futures, sizes, addrs, NULL);

    //bottom compartment can now be updated with children's contribution
    for (int c = 0; c < branchesCount; c++)
    {
        d[n-1]   -= values[c].b;
        rhs[n-1] -= values[c].rhs;
        hpx_lco_delete(futures[c], HPX_NULL);
    }

    double q;
    for (int i = n - 1; i >= 1; --i)
    {
        q = a[i]/d[i];
        d[i-1]   -= q * b[i];
        rhs[i-1] -= q * rhs[i];
    }

    delete [] futures;
    delete [] addrs;
    delete [] sizes;
    delete [] values;

    //value to be decremented will be sent to parent branch (except soma)
    if (!local->isSoma)
    {
        BackTriangFuture futureData;
        q = a[0] / d[0];
        futureData.b   = q * b[0];
        futureData.rhs = q * rhs[0];
        neurox_hpx_unpin_continue(futureData);
    }
    neurox_hpx_unpin;
}


hpx_action_t HinesSolver::gaussianFwdSubstitution = 0;
int HinesSolver::gaussianFwdSubstitution_handler(const double * parentRHS_ptr, const size_t)
{
    neurox_hpx_pin(Branch);
    double *b   = local->b;
    double *d   = local->d;
    double *rhs = local->rhs;
    int n = local->n;

    if(local->isSoma)
    {
        rhs[0] /= d[0];
    }
    else for (int i=0; i<n; i++)
    {
        rhs[i] -= b[i] * (i==0 ? *parentRHS_ptr : rhs[i-1]);
        rhs[i] /= d[i];
    }

    double childrenRHS=rhs[n-1];
    neurox_hpx_recursive_branch_sync(HinesSolver::gaussianFwdSubstitution, &childrenRHS, sizeof(childrenRHS));
    neurox_hpx_unpin;
}

hpx_action_t HinesSolver::setupMatrixLHS = 0;
int HinesSolver::setupMatrixLHS_handler()
{
    neurox_hpx_pin(Branch);
    int n = local->n;
    double *a   = local->a;
    double *b   = local->b;
    double *d   = local->d;
    int branchesCount = local->branchesCount;

    double returnValue = -1; //contribution to upper branch

    //for (i = i2; i < i3; ++i))
    if (!local->isSoma)
    {
        d[0] -= b[0];
        returnValue = a[0];
    }
    for (int i=1; i<local->n; i++)
    {
        d[i]   -= b[i];
        d[i-1] -= a[i];
    }

    //send/receive contribution to/from branches
    hpx_t * futures = new hpx_t[branchesCount];
    void  ** addrs  = new void*[branchesCount];
    size_t * sizes  = new size_t[branchesCount];
    double * values = new double[branchesCount];
    for (int c = 0; c<local->branchesCount; c++)
    {
        futures[c] = hpx_lco_future_new(sizeof(double));
        addrs[c]   = &values[c];
        sizes[c]   = sizeof(double);
        hpx_call(local->branches[c], HinesSolver::setupMatrixLHS, futures[c]);
    }

    if (branchesCount > 0) //required or fails
        hpx_lco_get_all(branchesCount, futures, sizes, addrs, NULL);

    //received contributions, can now update value
    for (int c = 0; c < branchesCount; c++)
    {
        d[n-1] -= values[c];
        hpx_lco_delete(futures[c], HPX_NULL);
    }

    delete [] futures;
    delete [] addrs;
    delete [] sizes;
    delete [] values;

    if (!local->isSoma)
        neurox_hpx_unpin_continue(returnValue);
    neurox_hpx_unpin;
}

void HinesSolver::registerHpxActions()
{
    neurox_hpx_register_action(1, HinesSolver::setupMatrixRHS);
    neurox_hpx_register_action(0, HinesSolver::setupMatrixLHS);
    neurox_hpx_register_action(0, HinesSolver::gaussianBackTriangulation);
    neurox_hpx_register_action(1, HinesSolver::gaussianFwdSubstitution);
}
