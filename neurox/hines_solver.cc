#include "neurox/neurox.h"

#include <algorithm>
#include <cstring>
#include <numeric>

using namespace neurox;
using namespace neurox::interpolators;

HinesSolver::~HinesSolver() {}

void HinesSolver::CommunicateConstants(const Branch *branch) {
  const floble_t *a = branch->nt_->_actual_a;
  const floble_t *b = branch->nt_->_actual_b;
  const int n = branch->nt_->end;
  Branch::BranchTree *branch_tree = branch->branch_tree_;

  if (!branch_tree) return;

  // all branches except top
  if (!branch->soma_) {
    floble_t to_parent_a_b[2] = {a[0], b[0]};
    floble_t from_parent_a_b[2];
    hpx_lco_set_rsync(branch_tree->with_parent_lco_[3], sizeof(floble_t) * 2,
                      &to_parent_a_b);
    hpx_lco_get_reset(branch_tree->with_parent_lco_[4], sizeof(floble_t) * 2,
                      &from_parent_a_b);
    branch_tree->parent_a_ = from_parent_a_b[0];
    branch_tree->parent_b_ = from_parent_a_b[1];
  }

  // all branches with leaves
  if (branch_tree != nullptr && branch_tree->branches_count_ > 0) {
    branch_tree->children_a_ = new floble_t[branch_tree->branches_count_];
    branch_tree->children_b_ = new floble_t[branch_tree->branches_count_];
    floble_t to_children_a_b[2] = {a[n - 1], b[n - 1]};
    floble_t from_child_a_b[2];
    for (offset_t c = 0; c < branch_tree->branches_count_; c++) {
      hpx_lco_get_reset(branch_tree->with_children_lcos_[c][3],
                        sizeof(floble_t) * 2, &from_child_a_b);
      branch_tree->children_a_[c] = from_child_a_b[0];
      branch_tree->children_b_[c] = from_child_a_b[1];
      hpx_lco_set_rsync(branch_tree->with_children_lcos_[c][4],
                        sizeof(floble_t) * 2, &to_children_a_b);
    }
  }
}

void HinesSolver::SynchronizeThresholdV(const Branch *branch,
                                        floble_t *threshold_v) {
  if (branch->soma_) {
    if (branch->thvar_ptr_)  // if I hold the value (Coreneuron base case)
      *threshold_v = *branch->thvar_ptr_;
    else
      // If not, wait for the value to be updated by AIS
      hpx_lco_get_reset(branch->branch_tree_->with_children_lcos_[0][0],
                        sizeof(floble_t), threshold_v);
  } else if (branch->thvar_ptr_)  // if AIS, send value to soma
    hpx_lco_set(branch->branch_tree_->with_parent_lco_[0], sizeof(floble_t),
                branch->thvar_ptr_, HPX_NULL, HPX_NULL);
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

  /* now the internal axial currents.
    The extracellular mechanism contribution is already done.
        rhs += ai_j*(vi_j - vi)
  */

  floble_t dv = 0;
  // first local future for downwards V, second for upwards A*dv
  if (!branch->soma_)  // all top compartments except soma
  {
    floble_t from_parent_v;  // get 'v[p[i]]' from parent
    hpx_lco_get_reset(branch_tree->with_parent_lco_[1], sizeof(floble_t),
                      &from_parent_v);

    dv = from_parent_v - v[0];
    rhs[0] -= b[0] * dv;

    floble_t to_parent_a = a[0] * dv;  // pass 'a[i]*dv' upwards to parent
    hpx_lco_set(branch_tree->with_parent_lco_[2], sizeof(floble_t),
                &to_parent_a, HPX_NULL, HPX_NULL);  // NOTE: async set
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
    floble_t to_children_v = v[n - 1];  // dv = v[p[i]] - v[i]
    for (offset_t c = 0; c < branch_tree->branches_count_; c++)
      hpx_lco_set(branch_tree->with_children_lcos_[c][1], sizeof(floble_t),
                  &to_children_v, HPX_NULL, HPX_NULL);
    // TODO: shall these async sets be freed later?!

    floble_t from_children_a;  // rhs[p[i]] += a[i]*dv
    for (offset_t c = 0; c < branch_tree->branches_count_; c++) {
      hpx_lco_get_reset(branch_tree->with_children_lcos_[c][2],
                        sizeof(floble_t), &from_children_a);
      rhs[n - 1] += from_children_a;
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

  /* bottom compartment, get D and RHS from children */
  if (branch_tree != nullptr && branch_tree->branches_count_ > 0) {
    const floble_t *children_a = branch_tree->children_a_;
    const floble_t *children_b = branch_tree->children_b_;
    floble_t from_child_d_rhs[2];
    for (offset_t c = 0; c < branch_tree->branches_count_; c++) {
      hpx_lco_get_reset(branch_tree->with_children_lcos_[c][3],
                        sizeof(floble_t) * 2, &from_child_d_rhs);
      //    d[n - 1] -= from_child_d_rhs[0];
      //    rhs[n - 1] -= from_child_d_rhs[1];

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
    // floble_t to_parent_d_rhs[2] = {a[0]/d[0]*b[0], a[0]/d[0]*rhs[0]};
    floble_t to_parent_d_rhs[2] = {d[0], rhs[0]};
    hpx_lco_set(branch_tree->with_parent_lco_[3], sizeof(floble_t) * 2,
                &to_parent_d_rhs, HPX_NULL, HPX_NULL);
  }
}

void HinesSolver::ForwardSubstitution(Branch *branch) {
  const offset_t n = branch->nt_->end;
  const floble_t *b = branch->nt_->_actual_b;
  const floble_t *d = branch->nt_->_actual_d;
  const offset_t *p = branch->nt_->_v_parent_index;
  floble_t *rhs = branch->nt_->_actual_rhs;
  const Branch::BranchTree *branch_tree = branch->branch_tree_;

  /* top compartment: get RHS from parent */
  if (!branch->soma_) {
    floble_t from_parent_d_rhs[2];
    hpx_lco_get_reset(branch_tree->with_parent_lco_[4], sizeof(floble_t)*2,
                      &from_parent_d_rhs);
    //floble_t & parent_d = from_parent_d_rhs[0];
    floble_t & parent_rhs = from_parent_d_rhs[1];
    rhs[0] -= b[0] * parent_rhs;
  }

  rhs[0] /= d[0];  //(Coreneuron base case)

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

  /* bottom compartment, send to children RHS[n-1] */
  if (branch_tree != nullptr) {
    floble_t to_children_d_rhs[2] = {d[n-1], rhs[n-1]};
    for (offset_t c = 0; c < branch_tree->branches_count_; c++)
      hpx_lco_set(branch_tree->with_children_lcos_[c][4], sizeof(floble_t)*2,
                  &to_children_d_rhs, HPX_NULL, HPX_NULL);
  }
}

void HinesSolver::ForwardTriangulationLinearCable(Branch *branch) {
  const offset_t n = branch->nt_->end;
  const floble_t *a = branch->nt_->_actual_a;
  const floble_t *b = branch->nt_->_actual_b;
  const offset_t *p = branch->nt_->_v_parent_index;
  floble_t *rhs = branch->nt_->_actual_rhs;
  floble_t *d = branch->nt_->_actual_d;
  const Branch::BranchTree *branch_tree = branch->branch_tree_;

  floble_t pp;
  assert(p);

  /* top compartment, get from parent D and RHS */
  if (!branch->soma_) {
    floble_t from_parent_d_rhs[2];
    hpx_lco_get_reset(branch_tree->with_parent_lco_[3], sizeof(floble_t) * 2,
                      &from_parent_d_rhs);

    floble_t &parent_d = from_parent_d_rhs[0];
    floble_t &parent_rhs = from_parent_d_rhs[1];

    pp = branch_tree->parent_b_ / parent_d;
    d[0] -= pp * branch_tree->parent_a_;
    rhs[0] -= pp * parent_rhs;
  }

  /* middle compartments, linear use case only */
  for (offset_t i = 0; i < n - 1; i++) {
    assert(i == p[i + 1]);  // make sure there's no branching
    pp = b[i] / d[i];
    d[i + 1] -= pp * a[i];
    rhs[i + 1] -= pp * rhs[i];
  }

  /* bottom compartment, send to children D and RHS */
  if (branch_tree != nullptr && branch_tree->branches_count_ > 0) {
    floble_t to_children_d_rhs[2] = {d[n - 1], rhs[n - 1]};
    for (offset_t c = 0; c < branch_tree->branches_count_; c++)
      hpx_lco_set(branch_tree->with_children_lcos_[c][3], sizeof(floble_t) * 2,
                  &to_children_d_rhs, HPX_NULL, HPX_NULL);
  }
}

void HinesSolver::BackwardSubstitutionLinearCable(Branch *branch) {
  const offset_t n = branch->nt_->end;
  const floble_t *a = branch->nt_->_actual_a;
  floble_t *rhs = branch->nt_->_actual_rhs;
  floble_t *d = branch->nt_->_actual_d;
  const Branch::BranchTree *branch_tree = branch->branch_tree_;

  /* bottom compartment: get D and RHS from children */
  if (branch_tree != nullptr) {
    floble_t from_children_d_rhs[2];
    for (offset_t c = 0; c < branch_tree->branches_count_; c++) {
      hpx_lco_get_reset(branch_tree->with_children_lcos_[c][4],
                        sizeof(floble_t) * 2, &from_children_d_rhs);
      floble_t &from_child_rhs = from_children_d_rhs[1];
      rhs[n - 1] -= a[n - 1] * from_child_rhs;
    }
  }

  rhs[n - 1] /= d[n - 1];  //(Coreneuron base case)

  /* middle compartments, linear use case only */
  for (offset_t i = n - 2; i >= 0; i--) {
    rhs[i] -= a[i] * rhs[i + 1];
    rhs[i] /= d[i];
  }

  /* top compartment, send to parent rhs[0] */
  if (!branch->soma_) {
    floble_t to_parent_rhs = rhs[0];
    hpx_lco_set(branch_tree->with_parent_lco_[4], sizeof(floble_t),
                &to_parent_rhs, HPX_NULL, HPX_NULL);
  }
}

void HinesSolver::BackwardTriangulationBottomSubsection(Branch *branch) {
  const offset_t n = branch->nt_->end;
  const floble_t *a = branch->nt_->_actual_a;
  const floble_t *b = branch->nt_->_actual_b;
  const offset_t *p = branch->nt_->_v_parent_index;
  floble_t *rhs = branch->nt_->_actual_rhs;
  floble_t *d = branch->nt_->_actual_d;
  const Branch::BranchTree *branch_tree = branch->branch_tree_;

  floble_t pp;
  assert(branch_tree != nullptr && branch_tree->branches_count_ == 0);
  assert(p);

  /* bottom and middle compartments */
  for (offset_t i = n - 1; i >= 1; i--) {
    pp = a[i] / d[i];
    d[p[i]] -= pp * b[i];
    rhs[p[i]] -= pp * rhs[i];
  }

  /* top compartment, send to parent D and RHS */
  if (!branch->soma_) {
    floble_t to_parent_d_rhs[2] = {d[0], rhs[0]};
    hpx_lco_set(branch_tree->with_parent_lco_[4], sizeof(floble_t) * 2,
                &to_parent_d_rhs, HPX_NULL, HPX_NULL);
  }
}

void HinesSolver::ForwardSubstitutionBottomSubsection(Branch *branch) {
  const offset_t n = branch->nt_->end;
  const floble_t *b = branch->nt_->_actual_b;
  const floble_t *d = branch->nt_->_actual_d;
  const offset_t *p = branch->nt_->_v_parent_index;
  floble_t *rhs = branch->nt_->_actual_rhs;
  const Branch::BranchTree *branch_tree = branch->branch_tree_;

  assert(branch_tree != nullptr && branch_tree->branches_count_ == 0);
  assert(p);

  /* top compartment: get D and RHS from parent */
  if (!branch->soma_) {
    floble_t from_parent_d_rhs[2];
    hpx_lco_get_reset(branch_tree->with_parent_lco_[3], sizeof(floble_t) * 2,
                      &from_parent_d_rhs);
    floble_t & parent_rhs = from_parent_d_rhs[1];
    rhs[0] -= b[0] * parent_rhs;
  }

  rhs[0] /= d[0];  //(Coreneuron base case)

  /* middle and bottom compartments */
  for (offset_t i = 1; i < n; i++) {
    rhs[i] -= b[i] * rhs[p[i]];  // reads from parent
    rhs[i] /= d[i];
  }
}

void HinesSolver::SolveTreeMatrix(Branch *branch) {
  /* Multisplit: neuron split into trees of linear cable if not a leaf of the
   * tree, or a branched subsection at the bottom of the tree otherwise.
   * Two-way Gaussian Elimination allows for parallelism of subsections. */
  if (branch->branch_tree_) {
    // if (branch->branch_tree_->branches_count_ > 0) {
    /* Regular Gaussian Elimination for tridiagonal matrix (linear cable) */
    //  HinesSolver::ForwardTriangulationLinearCable(branch);
    //  HinesSolver::BackwardSubstitutionLinearCable(branch);
    //} else {
    /* Inverted Gaussian Elimination for branched subsection */
    //  HinesSolver::BackwardTriangulationBottomSubsection(branch);
    //  HinesSolver::ForwardSubstitutionBottomSubsection(branch);
    //}
    // return;
  }

  /* Inverted Gaussian Elimination: from leaves to root
   * (or thread reduce+spwan version of branched Gaussian Elimination) */
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
