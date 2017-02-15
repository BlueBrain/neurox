#include "neurox/neurox.h"
#include <cstring>
#include <algorithm>
#include <numeric>

using namespace neurox;
using namespace neurox::Solver;

HinesSolver::~HinesSolver(){}

void HinesSolver::setupMatrixRHS(Branch * local)
{
    const offset_t n = local->nt->end;
    const floble_t *a   = local->nt->_actual_a;
    const floble_t *b   = local->nt->_actual_b;
    const floble_t *v   = local->nt->_actual_v;
    const offset_t *p = local->nt->_v_parent_index;
    floble_t *rhs = local->nt->_actual_rhs;
    const Branch::NeuronTree * neuronTree = local->neuronTree;

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
        //we use first local future for V, and second for RHS
        if (!local->soma) //all branches except top
        {
            floble_t fromParentV; //get 'v[p[i]]' from parent
            hpx_lco_get_reset(neuronTree->localLCO[0],
                    sizeof(floble_t), &fromParentV);
            dv = fromParentV - v[0];
            rhs[0] -= b[0]*dv;

            floble_t toParentRHS = a[0]*dv; //pass 'a[i]*dv' upwards to parent
            hpx_lco_set_rsync(neuronTree->localLCO[1],
                    sizeof(floble_t), &toParentRHS);
        }

        //middle compartments
        for (offset_t i=1; i<n; i++)
        {
            dv = v[i-1] - v[i];
            rhs[i] -= b[i]*dv;
            rhs[i-1] += a[i]*dv;
        }

        //bottom compartment
        floble_t toChildrenV = v[n-1];
        for (offset_t c=0; c<neuronTree->branchesCount; c++) //rhs[p[i]] += a[i]*dv
            hpx_lco_set_rsync(neuronTree->branchesLCOs[c][0],
                    sizeof(floble_t), &toChildrenV);

        floble_t fromChildrenRHS;
        for (offset_t c=0; c<neuronTree->branchesCount; c++) //rhs[p[i]] += a[i]*dv
        {
            hpx_lco_get_reset(neuronTree->branchesLCOs[c][1],
                    sizeof(floble_t), &fromChildrenRHS);
            rhs[n-1] += fromChildrenRHS;
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
    const Branch::NeuronTree * neuronTree = local->neuronTree;

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
    assert(0);
}

void HinesSolver::forwardSubstituion(Branch *local)
{
    assert(0);
}
