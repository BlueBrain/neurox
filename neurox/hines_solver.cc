#include "neurox/neurox.h"

#include <algorithm>
#include <cstring>
#include <numeric>

using namespace neurox;
using namespace neurox::interpolators;

HinesSolver::~HinesSolver() {}

// TODO: shall these async sets be freed later?!

void HinesSolver::InitializeBranchConstants(const Branch *branch) {
  const floble_t *a = branch->nt_->_actual_a;
  const floble_t *b = branch->nt_->_actual_b;
  Branch::BranchTree *bt = branch->branch_tree_;

  if (!bt) return;

  const int channel = 1;

  // all branches except top
  if (!branch->soma_) {
    floble_t to_parent_a_b[2] = {a[0], b[0]};
    assert(bt->with_children_lcos_[channel] != HPX_NULL);
    hpx_lco_set_rsync(bt->with_parent_lco_[channel], sizeof(floble_t) * 2,
                      &to_parent_a_b);
    bt->parent_v_ = input_params_->voltage_;
    bt->parent_rhs_ = -1;  // set by UpdateVoltagesWithRHS
  }

  // all branches with leaves
  if (bt != nullptr && bt->branches_count_ > 0) {
    bt->children_a_ = new floble_t[bt->branches_count_];
    bt->children_b_ = new floble_t[bt->branches_count_];
    bt->children_rhs_ = new floble_t[bt->branches_count_];
    bt->children_v_ = new floble_t[bt->branches_count_];
    floble_t from_child_a_b[2];
    for (offset_t c = 0; c < bt->branches_count_; c++) {
      assert(bt->with_children_lcos_[c][channel] != HPX_NULL);
      hpx_lco_get_reset(bt->with_children_lcos_[c][channel],
                        sizeof(floble_t) * 2, &from_child_a_b);
      bt->children_a_[c] = from_child_a_b[0];
      bt->children_b_[c] = from_child_a_b[1];
      bt->children_v_[c] = input_params_->voltage_;
      bt->children_rhs_[c] = -1;  // set by UpdateVoltagesWithRHS
    }
  }
}

void HinesSolver::SynchronizeThresholdV(const Branch *branch,
                                        floble_t *threshold_v) {
  const int channel = 0;
  const int ais_branch = 0;

  if (branch->soma_) {
    if (branch->thvar_ptr_)  // if I hold the value (Coreneuron base case)
      *threshold_v = *branch->thvar_ptr_;
    else
      // If not, wait for the value to be updated by AIS
      hpx_lco_get_reset(
          branch->branch_tree_->with_children_lcos_[ais_branch][channel],
          sizeof(floble_t), threshold_v);
  } else if (branch->thvar_ptr_)  // if AIS, send value to soma
    hpx_lco_set(branch->branch_tree_->with_parent_lco_[channel],
                sizeof(floble_t), branch->thvar_ptr_, HPX_NULL, HPX_NULL);
}

void HinesSolver::ResetArray(const Branch *branch, floble_t *arr) {
  std::fill(arr, arr + branch->nt_->end, 0.0);
}

void HinesSolver::SetupMatrixRHS(Branch *branch) {
  const offset_t n = branch->nt_->end;
  const floble_t *a = branch->nt_->_actual_a;
  const floble_t *b = branch->nt_->_actual_b;
  const floble_t *v = branch->nt_->_actual_v;
  const offset_t *p = branch->nt_->_v_parent_index;
  floble_t *rhs = branch->nt_->_actual_rhs;
  const Branch::BranchTree *branch_tree = branch->branch_tree_;
  const floble_t *children_a = branch_tree->children_a_;
  const floble_t *children_v = branch_tree->children_v_;

  /* now the internal axial currents.
    The extracellular mechanism contribution is already done.
        rhs += ai_j*(vi_j - vi)
  */

  floble_t dv = 0;
  // first local future for downwards V, second for upwards A*dv
  if (!branch->soma_)  // all top compartments except soma
  {
    dv = branch_tree->parent_v_ - v[0];
    rhs[0] -= b[0] * dv;
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
  if (branch_tree != nullptr && branch_tree->branches_count_ > 0) {
    for (offset_t c = 0; c < branch_tree->branches_count_; c++) {
      dv = v[n - 1] - children_v[c];
      rhs[n - 1] += children_a[c] * dv;
    }
  }
}

void HinesSolver::SetupMatrixDiagonal(Branch *branch) {
  const offset_t n = branch->nt_->end;
  const floble_t *a = branch->nt_->_actual_a;
  const floble_t *b = branch->nt_->_actual_b;
  const offset_t *p = branch->nt_->_v_parent_index;
  floble_t *d = branch->nt_->_actual_d;
  const Branch::BranchTree *branch_tree = branch->branch_tree_;

  // we use third local future for contribution A to parent's D
  if (!branch->soma_)  // all branches except top
  {
    d[0] -= b[0];
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
  if (branch_tree != nullptr && branch_tree->branches_count_ > 0) {
    for (offset_t c = 0; c < branch_tree->branches_count_; c++) {
      d[n - 1] -= branch_tree->children_a_[c];
    }
  }
}

void HinesSolver::BackwardTriangulation(Branch *branch) {
  const offset_t n = branch->nt_->end;
  const floble_t *a = branch->nt_->_actual_a;
  const floble_t *b = branch->nt_->_actual_b;
  const offset_t *p = branch->nt_->_v_parent_index;
  floble_t *rhs = branch->nt_->_actual_rhs;
  floble_t *d = branch->nt_->_actual_d;
  const Branch::BranchTree *branch_tree = branch->branch_tree_;

  floble_t pp;
  const int channel = 1;

  /* bottom compartment, get D and RHS from children */
  if (branch_tree != nullptr && branch_tree->branches_count_ > 0) {
    const floble_t *children_a = branch_tree->children_a_;
    const floble_t *children_b = branch_tree->children_b_;
    floble_t from_child_d_rhs[2];
    for (offset_t c = 0; c < branch_tree->branches_count_; c++) {
      hpx_lco_get_reset(branch_tree->with_children_lcos_[c][channel],
                        sizeof(floble_t) * 2, &from_child_d_rhs);

      floble_t &child_d = from_child_d_rhs[0];
      floble_t &child_rhs = from_child_d_rhs[1];
      pp = children_a[c] / child_d;
      d[n - 1] -= pp * children_b[c];
      rhs[n - 1] -= pp * child_rhs;
    }
  }

  /* middle compartments */
  if (p)  // may have bifurcations (Coreneuron base case)
    for (offset_t i = n - 1; i >= 1; i--) {
      pp = a[i] / d[i];
      d[p[i]] -= pp * b[i];
      rhs[p[i]] -= pp * rhs[i];
    }
  else  // a leaf on the tree
    for (offset_t i = n - 1; i >= 1; i--) {
      pp = a[i] / d[i];
      d[i - 1] -= pp * b[i];
      rhs[i - 1] -= pp * rhs[i];
    }

  /* top compartment, send to parent D and RHS */
  if (!branch->soma_) {
    floble_t to_parent_d_rhs[2] = {d[0], rhs[0]};
    hpx_lco_set(branch_tree->with_parent_lco_[channel], sizeof(floble_t) * 2,
                &to_parent_d_rhs, HPX_NULL, HPX_NULL);
  }
}

void HinesSolver::ForwardSubstitution(Branch *branch) {
  const offset_t n = branch->nt_->end;
  const floble_t *b = branch->nt_->_actual_b;
  const floble_t *d = branch->nt_->_actual_d;
  const offset_t *p = branch->nt_->_v_parent_index;
  floble_t *rhs = branch->nt_->_actual_rhs;
  Branch::BranchTree *branch_tree = branch->branch_tree_;

  const int channel_triang_rhs = 2;
  const int channel_final_rhs = 3;

  /* top compartment: */
  if (!branch->soma_) {
    /* get RHS from parent */
    floble_t &from_parent_rhs = branch_tree->parent_rhs_;
    hpx_lco_get_reset(branch_tree->with_parent_lco_[channel_triang_rhs],
                      sizeof(floble_t), &from_parent_rhs);

    rhs[0] -= b[0] * from_parent_rhs;
    rhs[0] /= d[0];

    /* pass final RHS to parent, will be needed by UpdateVoltagesWithRHS */
    floble_t &to_parent_rhs = rhs[0];
    hpx_lco_set(branch_tree->with_parent_lco_[channel_final_rhs],
                sizeof(floble_t), &to_parent_rhs, HPX_NULL, HPX_NULL);
  } else {
    rhs[0] /= d[0];  //(Coreneuron base case)
  }

  /* middle compartments */
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

  /* bottom compartment, send to children RHS[n-1] for their substitution,
     and get RHS[n-1] from children after they substituted */
  if (branch_tree != nullptr) {
    floble_t &to_children_rhs = rhs[n - 1];
    for (offset_t c = 0; c < branch_tree->branches_count_; c++)
      hpx_lco_set(branch_tree->with_children_lcos_[c][channel_triang_rhs],
                  sizeof(floble_t), &to_children_rhs, HPX_NULL, HPX_NULL);

    /* final RHS value from each children, after ForwardSubstitution */
    for (offset_t c = 0; c < branch_tree->branches_count_; c++)
      hpx_lco_get_reset(branch_tree->with_children_lcos_[c][channel_final_rhs],
                  sizeof(floble_t), &branch_tree->children_rhs_[c]);
  }
}

void HinesSolver::SolveTreeMatrix(Branch *branch) {
  /* Inverted Gaussian Elimination: from leaves to root */
  HinesSolver::BackwardTriangulation(branch);
  HinesSolver::ForwardSubstitution(branch);
}

void HinesSolver::UpdateVoltagesWithRHS(Branch *branch) {
  const floble_t *rhs = branch->nt_->_actual_rhs;
  floble_t *v = branch->nt_->_actual_v;

  // Reminder: after Gaussian Elimination, RHS is dV/dt (?)
  const floble_t second_order_multiplier = input_params_->second_order_ ? 2 : 1;
  for (int i = 0; i < branch->nt_->end; i++)
    v[i] += second_order_multiplier * rhs[i];

  // update branch-parallelism use case
  Branch::BranchTree *bt = branch->branch_tree_;
  if (bt != nullptr) {
    if (!branch->soma_)
      bt->parent_v_ += second_order_multiplier * bt->parent_rhs_;

    for (offset_t c = 0; c < bt->branches_count_; c++) {
      bt->children_v_[c] += second_order_multiplier * bt->children_rhs_[c];
    }
  }
}

void HinesSolver::ResetRHSandDNoCapacitors(Branch *branch) {
  const VariableTimeStep *vardt = (VariableTimeStep *)branch->interpolator_;

  floble_t *rhs = branch->nt_->_actual_rhs;
  floble_t *d = branch->nt_->_actual_d;
  const int *node_ids = vardt->no_cap_node_ids_;
  int nd = -1;
  for (int i = 0; i < vardt->no_cap_node_ids_count_; i++) {
    nd = node_ids[i];
    d[nd] = 0;
    rhs[nd] = 0;
  }
}

void HinesSolver::ResetRHSNoCapacitors(Branch *branch) {
  const VariableTimeStep *vardt = (VariableTimeStep *)branch->interpolator_;

  floble_t *rhs = branch->nt_->_actual_rhs;
  const int *node_ids = vardt->no_cap_node_ids_;
  for (int i = 0; i < vardt->no_cap_node_ids_count_; i++) {
    rhs[node_ids[i]] = 0;
  }
}

void HinesSolver::SetupMatrixVoltageNoCapacitors(Branch *branch) {
  const VariableTimeStep *vardt = (VariableTimeStep *)branch->interpolator_;

  const floble_t *a = branch->nt_->_actual_a;
  const floble_t *b = branch->nt_->_actual_b;
  const int *p = branch->nt_->_v_parent_index;
  const int *no_cap_child = vardt->no_cap_child_ids_;
  const int *no_cap_node = vardt->no_cap_node_ids_;

  floble_t *rhs = branch->nt_->_actual_rhs;
  floble_t *d = branch->nt_->_actual_d;
  floble_t *v = branch->nt_->_actual_v;
  int nd = -1, pnd = -1;

  // parent axial current
  for (int i = 0; i < vardt->no_cap_node_ids_count_; i++) {
    nd = no_cap_node[i];
    rhs[nd] += d[nd] * v[nd];
    if (nd > 0)  // has parent
    {
      rhs[nd] -= b[nd] * v[p[nd]];
      d[nd] -= b[nd];
    }
  }

  // child axial current (following from global v_parent)
  for (int i = 0; i < vardt->no_cap_child_ids_count_; i++) {
    nd = no_cap_child[i];
    pnd = p[nd];
    rhs[pnd] -= a[nd] * v[nd];
    d[pnd] -= a[nd];
  }

  for (int i = 0; i < vardt->no_cap_node_ids_count_; i++) {
    nd = no_cap_node[i];
    v[nd] = rhs[nd] / d[nd];
  }
  // no_cap voltages are now consistent with adjacent voltages
}
