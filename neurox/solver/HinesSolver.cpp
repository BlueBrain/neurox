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
        //first local future for downwards V, second for upwards A*dv
        if (!local->soma) //all top compartments except soma
        {
            floble_t fromParentV; //get 'v[p[i]]' from parent
            hpx_lco_get_reset(neuronTree->withParentLCO[0],
                    sizeof(floble_t), &fromParentV);

            dv = fromParentV - v[0];
            rhs[0] -= b[0]*dv;

            floble_t toParentA = a[0]*dv; //pass 'a[i]*dv' upwards to parent
            hpx_lco_set_rsync(neuronTree->withParentLCO[1],
                    sizeof(floble_t), &toParentA);
        }

        //middle compartments
        if (p) //may have bifurcations (Coreneuron base case)
          for (offset_t i=1; i<n; i++)
          {
            dv = v[p[i]] - v[i];  //reads from parent
            rhs[i] -= b[i]*dv;
            rhs[p[i]] += a[i]*dv; //writes to parent
          }
        else //a leaf on the tree
          for (offset_t i=1; i<n; i++)
          {
            dv = v[i-1] - v[i];
            rhs[i] -= b[i]*dv;
            rhs[i-1] += a[i]*dv;
          }

        //bottom compartment (when there is branching)
        if (neuronTree!=nullptr)
        {
          floble_t toChildrenV = v[n-1]; //dv = v[p[i]] - v[i]
          for (offset_t c=0; c<neuronTree->branchesCount; c++)
              hpx_lco_set_rsync(neuronTree->withChildrenLCOs[c][0],
                    sizeof(floble_t), &toChildrenV);

          floble_t fromChildrenA; //rhs[p[i]] += a[i]*dv
          for (offset_t c=0; c<neuronTree->branchesCount; c++)
          {
              hpx_lco_get_reset(neuronTree->withChildrenLCOs[c][1],
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

        //we use third local future for contribution A to parent's D
        if (!local->soma) //all branches except top
        {
            floble_t toParentA=-1;
            d[0] -= b[0];
            toParentA = a[0]; //pass 'a[i]' upwards to parent
            hpx_lco_set_rsync(neuronTree->withParentLCO[2],
                    sizeof(floble_t), &toParentA);
        }

        //middle compartments
        if (p) //may have bifurcations (Coreneuron base case)
          for (offset_t i=1; i<n; i++)
          {
            d[i] -= b[i];
            d[p[i]] -= a[i];
          }
        else //a leaf on the tree
          for (offset_t i=1; i<n; i++)
          {
            d[i] -= b[i];
            d[i-1] -= a[i];
          }

        //bottom compartment (when there is branching)
        if (neuronTree!=nullptr)
        {
          floble_t fromChildrenA;
          for (offset_t c=0; c<neuronTree->branchesCount; c++) //d[p[i]] -= a[i]
          {
            hpx_lco_get_reset(neuronTree->withChildrenLCOs[c][2],
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

        //bottom compartment (when there is branching)
        if (neuronTree!=nullptr)
        {
          floble_t fromChildrenD=-1, fromChildrenRHS=-1; //pp*b[i] and pp*rhs[i] from children
          for (offset_t c=0; c<neuronTree->branchesCount; c++)
          {
            hpx_lco_get_reset(neuronTree->withChildrenLCOs[c][0],
                    sizeof(floble_t), &fromChildrenD);
            hpx_lco_get_reset(neuronTree->withChildrenLCOs[c][1],
                    sizeof(floble_t), &fromChildrenRHS);
            d[n-1]   -= fromChildrenD;
            rhs[n-1] -= fromChildrenRHS;
          }
        }

        //middle compartments
        if (p) //may have bifurcations (Coreneuron base case)
            for (offset_t i=n-1; i>=1; i--)
            {
                pp = a[i] / d[i];
                d[p[i]]   -= pp * b[i];   //writes to parent
                rhs[p[i]] -= pp * rhs[i]; //writes to parent
            }
        else //a leaf on the tree
            for (offset_t i=n-1; i>=1; i--)
            {
                pp = a[i] / d[i];
                d[i-1]   -= pp * b[i];
                rhs[i-1] -= pp * rhs[i];
            }


        //top compartment
        //we use 1st and 2nd futures for contribution pp*b[i] and pp*rhs[i] to parent D and RHS
        if (!local->soma) //all branches except top
        {
          floble_t toParentD=-1, toParentRHS=-1;
          if (!local->soma) //all branches except top
          {
            pp = a[0]/d[0];
            toParentD   = pp*b[0];   //pass 'pp*b[i]' upwards to parent
            toParentRHS = pp*rhs[0]; //pass 'pp*rhs[i]' upwards to parent
            hpx_lco_set_rsync(neuronTree->withParentLCO[0],
                    sizeof(floble_t), &toParentD);
            hpx_lco_set_rsync(neuronTree->withParentLCO[1],
                    sizeof(floble_t), &toParentRHS);
          }
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

        //we use third local future for contribution RHS from parent
        if (!local->soma) //all branches except top
        {
            floble_t fromParentRHS; //get 'rhs[p[i]]' from parent
            hpx_lco_get_reset(neuronTree->withParentLCO[2],
                    sizeof(floble_t), &fromParentRHS);

            rhs[0] -= b[0] * fromParentRHS;
            rhs[0] /= d[0];
        }
        else //top compartment only
        {
            rhs[0] /= d[0];
        }

        //middle compartments
        if (p) //may have bifurcations (Coreneuron base case)
          for (offset_t i=1; i<n; i++)
          {
              rhs[i] -= b[i] * rhs[p[i]]; //reads from parent
              rhs[i] /= d[i];
          }
        else //a leaf on the tree
          for (offset_t i=1; i<n; i++)
          {
              rhs[i] -= b[i] * rhs[i-1];
              rhs[i] /= d[i];
          }

        //bottom compartment (when there is branching)
        if (neuronTree!=nullptr)
        {
          floble_t toChildrenRHS = rhs[n-1];   //rhs[i] -= b[i] * rhs[p[i]];
          for (offset_t c=0; c<neuronTree->branchesCount; c++)
              hpx_lco_set_rsync(neuronTree->withChildrenLCOs[c][2],
                    sizeof(floble_t), &toChildrenRHS);
        }
}
