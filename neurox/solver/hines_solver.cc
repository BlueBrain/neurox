#include "neurox/neurox.h"

#include <algorithm>
#include <cstring>
#include <numeric>

using namespace neurox;
using namespace neurox::solver;

HinesSolver::~HinesSolver() {}

void HinesSolver::SynchronizeThresholdV(Branch *branch, floble_t *threshold_v) {
  if (branch->soma_) {
    if (branch->thvar_ptr_)  // if I hold the value (Coreneuron base case)
      *threshold_v = *branch->thvar_ptr_;
    else
      // If not, wait for the value to be updated by AIS
      hpx_lco_get_reset(branch->branch_tree_->with_children_lcos_[0][6],
                        sizeof(floble_t), threshold_v);
  } else if (branch->thvar_ptr_)  // if AIS, send value to soma
    hpx_lco_set_rsync(branch->branch_tree_->with_parent_lco_[6],
                      sizeof(floble_t), branch->thvar_ptr_);
}

void HinesSolver::ResetMatrixRHSandD(Branch *branch) {
  floble_t *rhs = branch->nt_->_actual_rhs;
  floble_t *d = branch->nt_->_actual_d;
  const int n = branch->nt_->end;

  for (int i = 0; i < n; i++) {
    rhs[i] = 0;
    d[i] = 0;
  }
}

void HinesSolver::ResetMatrixV(Branch *branch) {
  floble_t *v = branch->nt_->_actual_v;
  const int n = branch->nt_->end;
  for (int i = 0; i < n; i++)
    v[i] = 0;
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
    hpx_lco_get_reset(branch_tree->with_parent_lco_[0], sizeof(floble_t),
                      &from_parent_v);

    dv = from_parent_v - v[0];
    rhs[0] -= b[0] * dv;

    floble_t to_parent_a = a[0] * dv;  // pass 'a[i]*dv' upwards to parent
    hpx_lco_set_rsync(branch_tree->with_parent_lco_[1], sizeof(floble_t),
                      &to_parent_a);
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
      hpx_lco_set_rsync(branch_tree->with_children_lcos_[c][0],
                        sizeof(floble_t), &to_children_v);

    floble_t from_children_a;  // rhs[p[i]] += a[i]*dv
    for (offset_t c = 0; c < branch_tree->branches_count_; c++) {
      hpx_lco_get_reset(branch_tree->with_children_lcos_[c][1],
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
    floble_t to_parent_a = a[0];  // pass 'a[i]' upwards to parent
    hpx_lco_set_rsync(branch_tree->with_parent_lco_[2], sizeof(floble_t),
                      &to_parent_a);
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
    floble_t from_children_a;
    for (offset_t c = 0; c < branch_tree->branches_count_;
         c++)  // d[p[i]] -= a[i]
    {
      hpx_lco_get_reset(branch_tree->with_children_lcos_[c][2],
                        sizeof(floble_t), &from_children_a);
      d[n - 1] -= from_children_a;
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

  // bottom compartment (when there is branching)
  if (branch_tree != nullptr && branch_tree->branches_count_ > 0) {
    floble_t from_children_b,
        from_children_rhs;  // pp*b[i] and pp*rhs[i] from children
    for (offset_t c = 0; c < branch_tree->branches_count_; c++) {
      hpx_lco_get_reset(branch_tree->with_children_lcos_[c][3],
                        sizeof(floble_t), &from_children_b);
      hpx_lco_get_reset(branch_tree->with_children_lcos_[c][4],
                        sizeof(floble_t), &from_children_rhs);
      d[n - 1] -= from_children_b;
      rhs[n - 1] -= from_children_rhs;
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
  if (!branch->soma_)  // all branches except top
  {
    floble_t to_parent_b = -1, to_parent_rhs = -1;
    if (!branch->soma_)  // all branches except top
    {
      pp = a[0] / d[0];
      to_parent_b = pp * b[0];      // pass 'pp*b[i]' upwards to parent
      to_parent_rhs = pp * rhs[0];  // pass 'pp*rhs[i]' upwards to parent
      hpx_lco_set_rsync(branch_tree->with_parent_lco_[3], sizeof(floble_t),
                        &to_parent_b);
      hpx_lco_set_rsync(branch_tree->with_parent_lco_[4], sizeof(floble_t),
                        &to_parent_rhs);
    }
  }
}

void HinesSolver::ForwardSubstituion(Branch *branch) {
  const offset_t n = branch->nt_->end;
  const floble_t *b = branch->nt_->_actual_b;
  const floble_t *d = branch->nt_->_actual_d;
  const offset_t *p = branch->nt_->_v_parent_index;
  floble_t *rhs = branch->nt_->_actual_rhs;
  const Branch::BranchTree *branch_tree = branch->branch_tree_;

  // we use third local future for contribution RHS from parent
  if (!branch->soma_)  // all branches except top
  {
    floble_t from_parent_rhs;  // get 'rhs[p[i]]' from parent
    hpx_lco_get_reset(branch_tree->with_parent_lco_[5], sizeof(floble_t),
                      &from_parent_rhs);

    rhs[0] -= b[0] * from_parent_rhs;
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
  if (branch_tree != nullptr) {
    floble_t to_children_rhs = rhs[n - 1];  // rhs[i] -= b[i] * rhs[p[i]];
    for (offset_t c = 0; c < branch_tree->branches_count_; c++)
      hpx_lco_set_rsync(branch_tree->with_children_lcos_[c][5],
                        sizeof(floble_t), &to_children_rhs);
  }
}

void HinesSolver::UpdateVoltagesWithRHS(Branch *branch) {
  const floble_t *rhs = branch->nt_->_actual_rhs;
  floble_t *v = branch->nt_->_actual_v;

  //Reminder: after Gaussian Elimination, RHS is dV/dt (?)
  const floble_t second_order_multiplier = input_params_->second_order_ ? 2 : 1;
  for (int i = 0; i < branch->nt_->end; i++)
    v[i] += second_order_multiplier * rhs[i];
}

void HinesSolver::ResetNoCapacitanceRHSandD(
        Branch *branch, void* branch_cvodes_ptr)
{
    interpolators::VariableTimeStep * branch_cvodes =
            (interpolators::VariableTimeStep*) branch_cvodes_ptr;

    floble_t *rhs = branch->nt_->_actual_rhs;
    floble_t *d = branch->nt_->_actual_d;
    for (int i=0; i<branch_cvodes->no_cap_count; i++)
    {
        int nd = branch_cvodes->no_cap_node[i];
        d[nd]=0;
        rhs[nd]=0;
    }
}


void HinesSolver::NoCapacitanceVoltage(
        Branch * branch, void * branch_cvodes_ptr)
{
    interpolators::VariableTimeStep * branch_cvodes =
            (interpolators::VariableTimeStep*) branch_cvodes_ptr;

    floble_t *rhs = branch->nt_->_actual_rhs;
    floble_t *d = branch->nt_->_actual_d;
    const floble_t *a = branch->nt_->_actual_a;
    const floble_t *b = branch->nt_->_actual_b;
    const int * p = branch->nt_->_v_parent_index;
    floble_t *v = branch->nt_->_actual_v;
    int nd=-1, pnd=-1;
    int * no_cap_child = branch_cvodes->no_cap_child;
    int * no_cap_node = branch_cvodes->no_cap_node;

    //parent axial current
    for (int i=0; i<branch_cvodes->no_cap_count; i++)
    {
        nd = no_cap_node[i];
        rhs[nd] += d[nd] * v[nd];
        if (nd>0) //has parent
        {
            rhs[nd] -= b[nd]*v[p[nd]];
            d[nd] -= b[nd];
        }
    }

    //child axial current (following from global v_parent)
    for (int i=0; i<branch_cvodes->no_cap_child_count; i++)
    {
        nd = no_cap_child[i];
        pnd = p[nd];
        rhs[pnd] -= a[nd]*v[nd];
        d[pnd] -= a[nd];
    }

    for (int i=0; branch_cvodes->no_cap_count; i++)
    {
        nd = no_cap_node[i];
        v[nd] = rhs[nd] / d[nd];
    }
    // no_cap v's are now consistent with adjacent v's
}
