#include "neurox/neurox.h"
#include <cstring>
#include <algorithm>
#include <numeric>

using namespace neurox;
using namespace neurox::Solver;

HinesSolver::~HinesSolver(){}

void HinesSolver::resetMatrixRHSandD(Branch * local)
{
    floble_t *rhs = local->nt->_actual_rhs;
    floble_t *d = local->nt->_actual_d;

    for (int i=0; i<local->nt->end; i++)
    {
        rhs[i]=0;
        d[i]=0;
    }
}

void HinesSolver::setupMatrixRHS(Branch * local)
{
    const offset_t n = local->nt->end;
    const floble_t *a   = local->nt->_actual_a;
    const floble_t *b   = local->nt->_actual_b;
    const floble_t *v   = local->nt->_actual_v;
    const offset_t *p = local->nt->_v_parent_index;
    floble_t *rhs = local->nt->_actual_rhs;
    const Branch::BranchTree * neuronTree = local->branchTree;

    /* now the internal axial currents.
      The extracellular mechanism contribution is already done.
          rhs += ai_j*(vi_j - vi)
    */

    floble_t dv=0;
    if (neuronTree==NULL)
    {
      assert(local->nt->_v_parent_index);
      for (offset_t i=1; i<n; i++)
      {
          dv = v[p[i]] - v[i];  //reads from parent
          rhs[i] -= b[i]*dv;
          rhs[p[i]] += a[i]*dv; //writes to parent
      }
    }
    else
    {
        //first local future for downwards V, second for upwards A*dv
        if (!local->soma) //all top compartments except soma
        {
            floble_t fromParentV; //get 'v[p[i]]' from parent
            hpx_lco_get_reset(neuronTree->localLCO[0],
                    sizeof(floble_t), &fromParentV);

            dv = fromParentV - v[0];
            rhs[0] -= b[0]*dv;

            floble_t toParentA = a[0]*dv; //pass 'a[i]*dv' upwards to parent
            hpx_lco_set_rsync(neuronTree->localLCO[1],
                    sizeof(floble_t), &toParentA);
        }

        //middle compartments
        for (offset_t i=1; i<n; i++)
        {
            dv = v[i-1] - v[i];
            rhs[i] -= b[i]*dv;
            rhs[i-1] += a[i]*dv;
        }

        //bottom compartment
        floble_t toChildrenV = v[n-1]; //dv = v[p[i]] - v[i]
        for (offset_t c=0; c<neuronTree->branchesCount; c++)
            hpx_lco_set_rsync(neuronTree->branchesLCOs[c][0],
                    sizeof(floble_t), &toChildrenV);

        floble_t fromChildrenA; //rhs[p[i]] += a[i]*dv
        for (offset_t c=0; c<neuronTree->branchesCount; c++)
        {
            hpx_lco_get_reset(neuronTree->branchesLCOs[c][1],
                    sizeof(floble_t), &fromChildrenA);
            rhs[n-1] += fromChildrenA;
        }
    }
}

void HinesSolver::setupMatrixDiagonal(Branch * local)
{
    const offset_t n = local->nt->end;
    const floble_t *a   = local->nt->_actual_a;
    const floble_t *b   = local->nt->_actual_b;
    const offset_t *p = local->nt->_v_parent_index;
    floble_t *d = local->nt->_actual_d;
    const Branch::BranchTree * neuronTree = local->branchTree;

    if (neuronTree==NULL)
    {
        assert(p);
        for (offset_t i=1; i<n; i++)
        {
            d[i] -= b[i];
            d[p[i]] -= a[i]; //writes to parent
        }
    }
    else
    {
        //we use third local future for A contribution to prent
        floble_t toParentA=-1;
        if (!local->soma) //all branches except top
        {
            d[0] -= b[0];
            toParentA = a[0]; //pass 'a[i]' upwards to parent
            hpx_lco_set_rsync(neuronTree->localLCO[2],
                    sizeof(floble_t), &toParentA);
        }

        //middle compartments
        for (offset_t i=1; i<n; i++)
        {
            d[i] -= b[i];
            d[i-1] -= a[i];
        }

        //bottom compartment
        floble_t fromChildrenA;
        for (offset_t c=0; c<neuronTree->branchesCount; c++) //d[p[i]] -= a[i]
        {
            hpx_lco_get_reset(neuronTree->branchesLCOs[c][2],
                    sizeof(floble_t), &fromChildrenA);
            d[n-1] -= fromChildrenA;
        }
    }
}

void HinesSolver::backwardTriangulation(Branch *local)
{
    const offset_t n = local->nt->end;
    const floble_t *a   = local->nt->_actual_a;
    const floble_t *b   = local->nt->_actual_b;
    const offset_t *p = local->nt->_v_parent_index;
    floble_t *rhs = local->nt->_actual_rhs;
    floble_t *d = local->nt->_actual_d;
    const Branch::BranchTree * neuronTree = local->branchTree;

    floble_t pp=0;
    if (neuronTree==NULL)
    {
      assert(local->nt->_v_parent_index);
      for (offset_t i=n-1; i>=1; i--)
      {
          pp = a[i] / d[i];
          d[p[i]] -= pp * b[i];     //write to parent
          rhs[p[i]] -= pp * rhs[i]; //write to parent
      }
    }
    else
    {
        assert(0);
    }
}

void HinesSolver::forwardSubstituion(Branch *local)
{
    const offset_t n = local->nt->end;
    const floble_t *b   = local->nt->_actual_b;
    const floble_t *d   = local->nt->_actual_d;
    const offset_t *p = local->nt->_v_parent_index;
    floble_t *rhs = local->nt->_actual_rhs;
    const Branch::BranchTree * neuronTree = local->branchTree;

    if (neuronTree==NULL)
    {
      assert(local->nt->_v_parent_index);
      rhs[0] /= d[0];

      for (offset_t i=1; i<n; i++)
      {
          rhs[i] -= b[i] * rhs[p[i]]; //reads from parent
          rhs[i] /= d[i];
      }
    }
    else
    {
        assert(0);
    }
}
