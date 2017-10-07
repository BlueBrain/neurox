#include "neurox/neurox.h"

#include <iostream>

using namespace std;
using namespace neurox;

void tools::Vectorizer::ConvertToSOA(Branch* b) {
  // NOTE: arrays with memory-aligned allocation:
  // ml->pdata, nt->data, ml->nodeindices, v->parent_index, and ml->shadow_*;
  // data needs to add gaps, pdata needs new offset values

  assert(memb_func);
  assert(LAYOUT == 0);

  // get total counts
  int N = b->nt_->end;
  size_t old_data_size = 6 * N;
  size_t new_data_size = 6 * SizeOf(N);
  assert(new_data_size % Vectorizer::kSOAPadding == 0);

  for (int m = 0; m < mechanisms_count_; m++) {
    b->mechs_instances_[m]._nodecount_padded =
        SizeOf(b->mechs_instances_[m].nodecount);
    new_data_size +=
        b->mechs_instances_[m]._nodecount_padded * mechanisms_[m]->data_size_;
    old_data_size +=
        b->mechs_instances_[m].nodecount * mechanisms_[m]->data_size_;
  }

  // old-to-new offset map
  std::vector<int> data_offsets(old_data_size, -99999);

  // old thvar offset
  int thvar_idx = b->thvar_ptr_ ? b->thvar_ptr_ - &b->nt_->_actual_v[0]
                                : -1;  // compartment id (typically 2)
  assert(thvar_idx == -1 || (thvar_idx >= 0 && thvar_idx < N));

  double* data_new = New<double>(new_data_size);
  double* data_old = b->nt_->_data;

  // add padding to data for RHS, D, A, B, V and area
  for (int i = 0; i < 6; i++)
    for (size_t j = 0; j < N; j++) {
      int old_offset = N * i + j;
      int new_offset = SizeOf(N) * i + j;
      data_new[new_offset] = b->nt_->_data[old_offset];
      data_offsets.at(old_offset) = new_offset;
    }

  b->nt_->_actual_rhs = &data_new[SizeOf(N) * 0];
  b->nt_->_actual_d = &data_new[SizeOf(N) * 1];
  b->nt_->_actual_a = &data_new[SizeOf(N) * 2];
  b->nt_->_actual_b = &data_new[SizeOf(N) * 3];
  b->nt_->_actual_v = &data_new[SizeOf(N) * 4];
  b->nt_->_actual_area = &data_new[SizeOf(N) * 5];

  // convert AP-threshold pointer
  b->thvar_ptr_ = b->thvar_ptr_ ? &b->nt_->_actual_v[thvar_idx] : nullptr;

  // add padding and converting AoS->SoA in nt->data and update ml->data for
  // mechs instances
  unsigned old_offset_acc = N * 6;
  unsigned new_offset_acc = SizeOf(N) * 6;

  for (int m = 0; m < neurox::mechanisms_count_; m++) {
    Memb_list* instances = &b->mechs_instances_[m];
    double* instance_data_new = &data_new[new_offset_acc];

    int total_pdata_size =
        instances->_nodecount_padded * mechanisms_[m]->pdata_size_;
    int* pdata_new = New<int>(total_pdata_size);
    int* pdata_old = instances->pdata;

    for (int n = 0; n < instances->nodecount; n++)  // for every node
    {
      for (size_t i = 0; i < mechanisms_[m]->data_size_;
           i++)  // for every variable
      {
        int old_offset = mechanisms_[m]->data_size_ * n + i;
        int new_offset = SizeOf(instances->nodecount) * i + n;
        data_new[new_offset_acc + new_offset] = instances->data[old_offset];
        assert(b->nt_->_data[old_offset_acc + old_offset] ==
               instances->data[old_offset]);
        data_offsets.at(old_offset_acc + old_offset) =
            new_offset_acc + new_offset;
      }

      for (size_t i = 0; i < mechanisms_[m]->pdata_size_;
           i++)  // for every pointer
      {
        // padding
        int old_offset = mechanisms_[m]->pdata_size_ * n + i;   // SoA
        int new_offset = instances->_nodecount_padded * i + n;  // AoS
        pdata_new[new_offset] = pdata_old[old_offset];

        // get correct pdata offset: without branching, offsets are already
        // correct for both LAYOUTs and padding
        if (input_params_->branch_parallelism_depth_ > 0) {
          int ptype = memb_func[mechanisms_[m]->type_].dparam_semantics[i];
          bool is_pointer = ptype == -1 || (ptype > 0 && ptype < 1000);
          if (is_pointer)  // true for pointer to area in nt->data, or ion
                           // instance data
          {
            pdata_new[new_offset] =
                data_offsets.at(pdata_old[old_offset]);  // point to new offset
            assert(data_new[pdata_new[new_offset]] ==
                   data_old[pdata_old[old_offset]]);
          }
        }
        assert(pdata_new[new_offset] != -99999);
      }
    }

    old_offset_acc += mechanisms_[m]->data_size_ * instances->nodecount;
    new_offset_acc += mechanisms_[m]->data_size_ * SizeOf(instances->nodecount);

    // all instances processed, replace pointer by new padded data
    Delete(instances->pdata);
    instances->pdata = pdata_new;
    instances->data = instance_data_new;
  }

  // convert VecPlay continuous pointers
  for (int v = 0; v < b->nt_->n_vecplay; v++) {
    VecplayContinuousX* vc = (VecplayContinuousX*)b->nt_->_vecplay[v];
    int pd_offset = vc->pd_ - &b->nt_->_data[0];
    assert(pd_offset >= 0 && pd_offset <= b->nt_->_ndata);
    int pd_offset_new = data_offsets.at(pd_offset);
    vc->pd_ = &data_new[pd_offset_new];
  }

  Delete(b->nt_->_data);
  b->nt_->_data = data_new;
  b->nt_->_ndata = new_data_size;
}

void tools::Vectorizer::CallVecFunction(cvode_f_t func, NrnThread* nt,
                                        Memb_list* ml, int type) {
  double* p = nullptr;
  Datum* ppvar = nullptr;
  double v = 0.0;
  int* ni = ml->nodeindices;
  int cntml_actual = ml->nodecount;
  int cntml_padded = ml->_nodecount_padded;
  ThreadDatum* thread = ml->_thread;
  double* vec_v = nt->_actual_v;
  int nd_idx = -1;

  const int psize = GetMechanismFromType(type)->data_size_;
  const int ppsize = GetMechanismFromType(type)->pdata_size_;

#if LAYOUT == 1 /*AoS*/
  for (int iml = 0; iml < cntml_actual; ++iml) {
    p = ml->data + iml * psize;
    ppvar = ml->pdata + iml * ppsize;
#elif LAYOUT == 0 /*SoA*/
  p = ml->data;
  ppvar = ml->pdata;
  for (int iml = 0; iml < cntml_actual; ++iml) {
#endif
    nd_idx = ni[iml];
    v = vec_v[nd_idx];
    (*func)(iml, cntml_padded, p, ppvar, thread, nt, v);
  }
}
