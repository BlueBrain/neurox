#include "neurox/neurox.h"

#include <algorithm>
#include <cstring>
#include <numeric>

using namespace neurox;
using namespace neurox::solver;

HinesSolver::~HinesSolver() {}

void HinesSolver::SynchronizeThresholdV(Branch *local, floble_t *thresholdV) {
  if (local->soma_) {
    if (local->thvar_ptr_)  // if I hold the value (Coreneuron base case)
      *thresholdV = *local->thvar_ptr_;
    else
      // If not, wait for the value to be updated by AIS
      hpx_lco_get_reset(local->branch_tree_->with_children_lcos_[0][6],
                        sizeof(floble_t), thresholdV);
  } else if (local->thvar_ptr_)  // if AIS, send value to soma
    hpx_lco_set_rsync(local->branch_tree_->with_parent_lco_[6], sizeof(floble_t),
                      local->thvar_ptr_);
}

void HinesSolver::ResetMatrixRHSandD(Branch *local) {
  floble_t *rhs = local->nt_->_actual_rhs;
  floble_t *d = local->nt_->_actual_d;

  for (int i = 0; i < local->nt_->end; i++) {
    rhs[i] = 0;
    d[i] = 0;
  }
}

void HinesSolver::SetupMatrixRHS(Branch *local) {
  const offset_t n = local->nt_->end;
  const floble_t *a = local->nt_->_actual_a;
  const floble_t *b = local->nt_->_actual_b;
  const floble_t *v = local->nt_->_actual_v;
  const offset_t *p = local->nt_->_v_parent_index;
  floble_t *rhs = local->nt_->_actual_rhs;
  const Branch::BranchTree *neuronTree = local->branch_tree_;

  /* now the internal axial currents.
    The extracellular mechanism contribution is already done.
        rhs += ai_j*(vi_j - vi)
  */

  floble_t dv = 0;
  // first local future for downwards V, second for upwards A*dv
  if (!local->soma_)  // all top compartments except soma
  {
    floble_t fromParentV;  // get 'v[p[i]]' from parent
    hpx_lco_get_reset(neuronTree->with_parent_lco_[0], sizeof(floble_t),
                      &fromParentV);

    dv = fromParentV - v[0];
    rhs[0] -= b[0] * dv;

    floble_t toParentA = a[0] * dv;  // pass 'a[i]*dv' upwards to parent
    hpx_lco_set_rsync(neuronTree->with_parent_lco_[1], sizeof(floble_t),
                      &toParentA);
  }

  // middle compartments
  if (p)  // may have bifurcations (Coreneuron base case)
    for (offset_t i = 1; i < n; i++) {
      dv = v[p[i]] - v[i];  // reads from parent
      rhs[i] -= b[i] * dv;
      rhs[p[i]] += a[i] * dv;  // writes to parent
    }
  else  // a leaf on the tree
    for (offset_t i = 1; i < n; i++) {
      dv = v[i - 1] - v[i];
      rhs[i] -= b[i] * dv;
      rhs[i - 1] += a[i] * dv;
    }

  // bottom compartment (when there is branching)
  if (neuronTree != nullptr && neuronTree->branches_count_ > 0) {
    floble_t toChildrenV = v[n - 1];  // dv = v[p[i]] - v[i]
    for (offset_t c = 0; c < neuronTree->branches_count_; c++)
      hpx_lco_set_rsync(neuronTree->with_children_lcos_[c][0], sizeof(floble_t),
                        &toChildrenV);

    floble_t fromChildrenA;  // rhs[p[i]] += a[i]*dv
    for (offset_t c = 0; c < neuronTree->branches_count_; c++) {
      hpx_lco_get_reset(neuronTree->with_children_lcos_[c][1], sizeof(floble_t),
                        &fromChildrenA);
      rhs[n - 1] += fromChildrenA;
    }
  }
}

void HinesSolver::SetupMatrixDiagonal(Branch *local) {
  const offset_t n = local->nt_->end;
  const floble_t *a = local->nt_->_actual_a;
  const floble_t *b = local->nt_->_actual_b;
  const offset_t *p = local->nt_->_v_parent_index;
  floble_t *d = local->nt_->_actual_d;
  const Branch::BranchTree *neuronTree = local->branch_tree_;

  // we use third local future for contribution A to parent's D
  if (!local->soma_)  // all branches except top
  {
    d[0] -= b[0];
    floble_t toParentA = a[0];  // pass 'a[i]' upwards to parent
    hpx_lco_set_rsync(neuronTree->with_parent_lco_[2], sizeof(floble_t),
                      &toParentA);
  }

  // middle compartments
  if (p)  // may have bifurcations (Coreneuron base case)
    for (offset_t i = 1; i < n; i++) {
      d[i] -= b[i];
      d[p[i]] -= a[i];
    }
  else  // a leaf on the tree
    for (offset_t i = 1; i < n; i++) {
      d[i] -= b[i];
      d[i - 1] -= a[i];
    }

  // bottom compartment (when there is branching)
  if (neuronTree != nullptr && neuronTree->branches_count_ > 0) {
    floble_t fromChildrenA;
    for (offset_t c = 0; c < neuronTree->branches_count_;
         c++)  // d[p[i]] -= a[i]
    {
      hpx_lco_get_reset(neuronTree->with_children_lcos_[c][2], sizeof(floble_t),
                        &fromChildrenA);
      d[n - 1] -= fromChildrenA;
    }
  }
}

void HinesSolver::BackwardTriangulation(Branch *local) {
  const offset_t n = local->nt_->end;
  const floble_t *a = local->nt_->_actual_a;
  const floble_t *b = local->nt_->_actual_b;
  const offset_t *p = local->nt_->_v_parent_index;
  floble_t *rhs = local->nt_->_actual_rhs;
  floble_t *d = local->nt_->_actual_d;
  const Branch::BranchTree *neuronTree = local->branch_tree_;

  floble_t pp;

  // bottom compartment (when there is branching)
  if (neuronTree != nullptr && neuronTree->branches_count_ > 0) {
    floble_t fromChildrenB,
        fromChildrenRHS;  // pp*b[i] and pp*rhs[i] from children
    for (offset_t c = 0; c < neuronTree->branches_count_; c++) {
      hpx_lco_get_reset(neuronTree->with_children_lcos_[c][3], sizeof(floble_t),
                        &fromChildrenB);
      hpx_lco_get_reset(neuronTree->with_children_lcos_[c][4], sizeof(floble_t),
                        &fromChildrenRHS);
      d[n - 1] -= fromChildrenB;
      rhs[n - 1] -= fromChildrenRHS;
    }
  }

  // middle compartments
  if (p)  // may have bifurcations (Coreneuron base case)
    for (offset_t i = n - 1; i >= 1; i--) {
      pp = a[i] / d[i];
      d[p[i]] -= pp * b[i];      // writes to parent
      rhs[p[i]] -= pp * rhs[i];  // writes to parent
    }
  else  // a leaf on the tree
    for (offset_t i = n - 1; i >= 1; i--) {
      pp = a[i] / d[i];
      d[i - 1] -= pp * b[i];
      rhs[i - 1] -= pp * rhs[i];
    }

  // top compartment
  // we use 1st and 2nd futures for contribution pp*b[i] and pp*rhs[i] to parent
  // D and RHS
  if (!local->soma_)  // all branches except top
  {
    floble_t toParentB = -1, toParentRHS = -1;
    if (!local->soma_)  // all branches except top
    {
      pp = a[0] / d[0];
      toParentB = pp * b[0];      // pass 'pp*b[i]' upwards to parent
      toParentRHS = pp * rhs[0];  // pass 'pp*rhs[i]' upwards to parent
      hpx_lco_set_rsync(neuronTree->with_parent_lco_[3], sizeof(floble_t),
                        &toParentB);
      hpx_lco_set_rsync(neuronTree->with_parent_lco_[4], sizeof(floble_t),
                        &toParentRHS);
    }
  }
}

void HinesSolver::ForwardSubstituion(Branch *local) {
  const offset_t n = local->nt_->end;
  const floble_t *b = local->nt_->_actual_b;
  const floble_t *d = local->nt_->_actual_d;
  const offset_t *p = local->nt_->_v_parent_index;
  floble_t *rhs = local->nt_->_actual_rhs;
  const Branch::BranchTree *neuronTree = local->branch_tree_;

  // we use third local future for contribution RHS from parent
  if (!local->soma_)  // all branches except top
  {
    floble_t fromParentRHS;  // get 'rhs[p[i]]' from parent
    hpx_lco_get_reset(neuronTree->with_parent_lco_[5], sizeof(floble_t),
                      &fromParentRHS);

    rhs[0] -= b[0] * fromParentRHS;
    rhs[0] /= d[0];
  } else  // top compartment only (Coreneuron base case)
  {
    rhs[0] /= d[0];
  }

  // middle compartments
  if (p)  // may have bifurcations (Coreneuron base case)
    for (offset_t i = 1; i < n; i++) {
      rhs[i] -= b[i] * rhs[p[i]];  // reads from parent
      rhs[i] /= d[i];
    }
  else  // a leaf on the tree
    for (offset_t i = 1; i < n; i++) {
      rhs[i] -= b[i] * rhs[i - 1];
      rhs[i] /= d[i];
    }

  // bottom compartment (when there is branching)
  if (neuronTree != nullptr) {
    floble_t toChildrenRHS = rhs[n - 1];  // rhs[i] -= b[i] * rhs[p[i]];
    for (offset_t c = 0; c < neuronTree->branches_count_; c++)
      hpx_lco_set_rsync(neuronTree->with_children_lcos_[c][5], sizeof(floble_t),
                        &toChildrenRHS);
  }
}

void HinesSolver::UpdateV(Branch *local) {
  floble_t *rhs = local->nt_->_actual_rhs;
  floble_t *v = local->nt_->_actual_v;

  floble_t secondOrderMultiplier = input_params_->second_order_ ? 2 : 1;
  for (int i = 0; i < local->nt_->end; i++)
    v[i] += secondOrderMultiplier * rhs[i];
}
