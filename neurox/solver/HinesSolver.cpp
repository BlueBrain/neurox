#include "neurox/neurox.h"
#include <cstring>
#include <algorithm>
#include <numeric>

using namespace neurox;
using namespace neurox::Solver;

hpx_action_t HinesSolver::setupMatrixRHS = 0;
int HinesSolver::setupMatrixRHS_handler()
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(HinesSolver::setupMatrixRHS);
    for (int i=0; i<local->n; i++)
    {
        local->rhs[i]=0;
        local->d[i]=0;
    }
    local->callModFunction2(Mechanism::ModFunction::before_breakpoint);
    local->callModFunction2(Mechanism::ModFunction::current);
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t HinesSolver::gaussianFwdTriangulation = 0;
int HinesSolver::gaussianFwdTriangulation_handler(const double * vFromParent_ptr, const size_t)
{
  neurox_hpx_pin(Branch);
  const int n = local->n;
  const double *a   = local->a;
  const double *b   = local->b;
  const double *d   = local->d;
  const double *v   = local->v;
  double *rhs = local->rhs;
  const int *p = local->p;
  const int branchesCount = local->branchesCount;

  /* now the internal axial currents.
    The extracellular mechanism contribution is already done.
        rhs += ai_j*(vi_j - vi)
  */

  if (local->p!=NULL)
  {
    double dv=0;
    for (int i=1; i<n; i++)
    {
        dv = v[p[i]] - v[i];  //reads from parent
        rhs[i] -= b[i]*dv;
        rhs[p[i]] += a[i]*dv; //writes to parent
    }
  }
  else //multiSpliX
  {
    hpx_t * futures = branchesCount ? new hpx_t[branchesCount]  : nullptr;
    void  ** addrs  = branchesCount ? new void*[branchesCount]  : nullptr;
    size_t * sizes  = branchesCount ? new size_t[branchesCount] : nullptr;
    double * childrenValues = branchesCount ? new double[branchesCount] : nullptr;

    double valueForParent=-1;
    for (int i = local->isSoma ? 1 : 0 ; i <n; i++)
    {   //reads V from parent:
        double dv = (i==0 ? *vFromParent_ptr : v[i-1]) - v[i];
        rhs[i] -= b[i] *dv;
        if (i>1)
            rhs[i-1] += a[i]*dv;
        else
            valueForParent = a[i]*dv; //value to pass to parent
    }

    for (int c = 0; c < branchesCount; c++)
    {   //passes V to children and waits for their values of a[i]*dv
        futures[c] = hpx_lco_future_new(sizeof (double));
        addrs[c]   = &childrenValues[c];
        sizes[c]   = sizeof(double);
        hpx_call(local->branches[c], HinesSolver::gaussianFwdTriangulation, futures[c], &valueForParent, sizeof(valueForParent));
    }

    if (branchesCount > 0) //required or fails
        hpx_lco_get_all(branchesCount, futures, sizes, addrs, NULL);

    //bottom compartment can now be updated with children's contribution
    for (int c = 0; c < branchesCount; c++)
    {
        rhs[n-1] += childrenValues[c];
        hpx_lco_delete(futures[c], HPX_NULL);
    }

    delete [] futures;
    delete [] addrs;
    delete [] sizes;
    delete [] childrenValues;

    if (!local->isSoma) //send parent its value
        neurox_hpx_unpin_continue(valueForParent);
  }
  neurox_hpx_unpin;
}


hpx_action_t HinesSolver::setupMatrixLHS = 0;
int HinesSolver::setupMatrixLHS_handler()
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(HinesSolver::setupMatrixLHS);
    // calculate left hand side of
    //cm*dvm/dt = -i(vm) + is(vi) + ai_j*(vi_j - vi)
    //cx*dvx/dt - cm*dvm/dt = -gx*(vx - ex) + i(vm) + ax_j*(vx_j - vx)
    //with a matrix so that the solution is of the form [dvm+dvx,dvx] on the right
    //hand side after solving.
    //This is a common operation for fixed step, cvode, and daspk methods
    local->callModFunction2(Mechanism::ModFunction::jacob);

    //finitialize.c:nrn_finitialize()->set_tree_matrix_minimal->nrn_rhs (treeset_core.c)
    //now the cap current can be computed because any change to cm
    //by another model has taken effect.
    local->callModFunction2(Mechanism::ModFunction::jacobCapacitance);

    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}


hpx_action_t HinesSolver::gaussianBackSubstitution = 0;
int HinesSolver::gaussianBackSubstitution_handler()
{
    neurox_hpx_pin(Branch);
    const double *b   = local->b;
    const double *a   = local->a;
    const double *rhs = local->rhs;
    const double *v = local->v;
    const int *p = local->p;
    const int n = local->n;
    const int branchesCount = local->branchesCount;
    double *d   = local->d;

    if (p!=NULL)
    {
        for (int i=1; i<n; i++)
        {
            d[i] -= b[i];
            d[p[i]] -= a[i];
        }
    }
    else //multiSpliX
    {
      hpx_t * futures = branchesCount ? new hpx_t[branchesCount]  : nullptr;
      void  ** addrs  = branchesCount ? new void*[branchesCount]  : nullptr;
      size_t * sizes  = branchesCount ? new size_t[branchesCount] : nullptr;
      double * childrenValues = branchesCount ? new double[branchesCount] : nullptr;

      double valueForParent=-1;
      for (int i = local->isSoma ? 1 : 0 ; i <n; i++)
      {
          d[i] -= b[i];
          if (i>1)
              d[i-1] -= a[i];
          else
              valueForParent = a[i]; //value to pass to parent
      }

      for (int c = 0; c < branchesCount; c++)
      {   //waits for their values of a[i]
          futures[c] = hpx_lco_future_new(sizeof (double));
          addrs[c]   = &childrenValues[c];
          sizes[c]   = sizeof(double);
          hpx_call(local->branches[c], HinesSolver::gaussianBackSubstitution, futures[c]);
      }

      if (branchesCount > 0) //required or fails
          hpx_lco_get_all(branchesCount, futures, sizes, addrs, NULL);

      //bottom compartment can now be updated with children's contribution
      for (int c = 0; c < branchesCount; c++)
      {
          d[n-1] -= childrenValues[c];
          hpx_lco_delete(futures[c], HPX_NULL);
      }

      delete [] futures;
      delete [] addrs;
      delete [] sizes;
      delete [] childrenValues;
    }
    neurox_hpx_unpin;
}

void HinesSolver::registerHpxActions()
{
    neurox_hpx_register_action(0, HinesSolver::setupMatrixRHS);
    neurox_hpx_register_action(0, HinesSolver::setupMatrixLHS);
    neurox_hpx_register_action(1, HinesSolver::gaussianFwdTriangulation);
    neurox_hpx_register_action(0, HinesSolver::gaussianBackSubstitution);
}
