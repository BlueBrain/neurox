#include "neurox/neurox.h"
#include <cstring>
#include <algorithm>
#include <numeric>

using namespace neurox;
using namespace neurox::Solver;

HinesSolver::~HinesSolver(){}

void HinesSolver::forwardTriangulation(Branch * local)
{
    const int n = local->n;
    const double *a   = local->a;
    const double *b   = local->b;
    const double *v   = local->v;
    const int *p = local->p;
    double *rhs = local->rhs;
    const Branch::NeuronTree * neuronTree = local->neuronTree;

    /* now the internal axial currents.
      The extracellular mechanism contribution is already done.
          rhs += ai_j*(vi_j - vi)
    */

    double dv=0;
    if (neuronTree==NULL)
    {
      assert(local->p);
      for (int i=1; i<n; i++)
      {
          dv = v[p[i]] - v[i];  //reads from parent
          rhs[i] -= b[i]*dv;
          rhs[p[i]] += a[i]*dv; //writes to parent
      }
    }
    else
    {
        //we use first local future for V, and second for RHS

        if (!local->soma) //all branches except top
        {
            double fromParentV; //get 'v[p[i]]' from parent
            hpx_lco_get_reset(neuronTree->localLCO[0], sizeof(double), &fromParentV);
            dv = fromParentV - v[0];
            rhs[0] -= b[0]*dv;

            double toParentRHS = a[0]*dv; //pass 'a[i]*dv' upwards to parent
            hpx_lco_set_rsync(neuronTree->localLCO[1], sizeof(double), &toParentRHS);
        }

        //middle compartments
        for (int i=1; i<n-1; i++)
        {
            dv = v[i-1] - v[i];
            rhs[i] -= b[i]*dv;
            rhs[i-1] += a[i]*dv;
        }

        //bottom compartment
        assert(n>=2); //will fail otherwise
        dv = v[n-2] - v[n-1];  //dv = v[p[i]] - v[i];
        rhs[n-1] -= b[n-1]*dv; //rhs[i] -= b[i]*dv

        double toChildrenV = v[n-1];
        for (int c=0; c<neuronTree->branchesCount; c++) //rhs[p[i]] += a[i]*dv
            hpx_lco_set_rsync(neuronTree->branchesLCOs[c*2+0], sizeof(double), &toChildrenV);

        double fromChildrenRHS;
        for (int c=0; c<neuronTree->branchesCount; c++) //rhs[p[i]] += a[i]*dv
        {
            hpx_lco_get_reset(neuronTree->branchesLCOs[c*2+1], sizeof(double), &fromChildrenRHS);
            rhs[n-1] += fromChildrenRHS;
        }
    }
}

void HinesSolver::backSubstitution(Branch * local)
{
    const int n = local->n;
    const double *a   = local->a;
    const double *b   = local->b;
    const double *v   = local->v;
    const double *rhs = local->rhs;
    const int *p = local->p;
    double *d   = local->d;
    const Branch::NeuronTree * neuronTree = local->neuronTree;

    if (neuronTree==NULL)
    {
        assert(local->p);
        for (int i=1; i<n; i++)
        {
            d[i] -= b[i];
            d[p[i]] -= a[i]; //writes to parent
        }
    }
    else
    {
        double toParentA=-1;
        if (!local->soma) //all branches except top
        {
            d[0] -= b[0];
            toParentA = a[0]; //pass 'a[i]' upwards to parent
            hpx_lco_set_rsync(neuronTree->localLCO[1], sizeof(double), &toParentA);
        }

        //middle compartments
        for (int i=1; i<n-1; i++)
        {
            d[i] -= b[i];
            d[i-1] -= a[i];
        }

        //bottom compartment
        d[n-1] -= b[n-1]; //d[i] -= b[i]
        double fromChildrenA;
        for (int c=0; c<neuronTree->branchesCount; c++) //d[p[i]] -= a[i]
        {
            hpx_lco_get_reset(neuronTree->branchesLCOs[c*2+1], sizeof(double), &fromChildrenA);
            d[n-1] -= fromChildrenA;
        }
    }
}
