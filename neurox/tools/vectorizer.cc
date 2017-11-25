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
        if (input_params_->branch_parallelism_ ) {
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

void tools::Vectorizer::GroupBranchInstancesByCapacitors(
    const Branch* branch,              // in
    Memb_list** ml_no_capacitors_ptr,  // out (optional)
    Memb_list** ml_capacitors_ptr,     // out (optional)
    std::set<int>* capacitor_ids_ptr   // in  (optional)
    ) {
  // if not provided, build list of capacitor ids
  std::set<int> new_capacitor_ids;
  if (capacitor_ids_ptr == nullptr) {
    Memb_list* capac_instances =
        &branch->mechs_instances_[mechanisms_map_[CAP]];
    // get list of all nodes that are capacitors
    for (int c = 0; c < capac_instances->nodecount; c++) {
      int compartment_id = capac_instances->nodeindices[c];
      new_capacitor_ids.insert(compartment_id);
    }
  }

  std::set<int>& capacitor_ids =
      capacitor_ids_ptr ? *capacitor_ids_ptr : new_capacitor_ids;

  /* occvode.cpp::new_no_cap_memb(): get Memb_list for non-capacitor
  nodes only: pointers will point to same place in nt->data, we
  will re-order Memb_list to have no-caps first, and then
  update nodecount for no-caps instance to cover no-caps only */
  Memb_list* ml_no_capacitors = new Memb_list[neurox::mechanisms_count_];
  memcpy(ml_no_capacitors, branch->mechs_instances_,
         neurox::mechanisms_count_ * sizeof(Memb_list));

  // ml_capacitors is optional
  if (ml_capacitors_ptr != nullptr) {
    Memb_list*& ml_capacitors = *ml_capacitors_ptr;
    ml_capacitors = new Memb_list[neurox::mechanisms_count_];
    memcpy(ml_capacitors, branch->mechs_instances_,
           neurox::mechanisms_count_ * sizeof(Memb_list));
  }

  int total_data_offset = Vectorizer::SizeOf(branch->nt_->end) * 6;
  map<int, map<int, int>> ions_data_map;

  for (int m = 0; m < neurox::mechanisms_count_; m++) {
    Mechanism* mech = neurox::mechanisms_[m];
    Memb_list* instances = &branch->mechs_instances_[m];

    // neuron: "only point processes with currents are possibilities"
    bool mech_valid_in_phase_1 = mech->pnt_map_ && mech->memb_func_.current;

    int n_new = 0;
    int data_size = mech->data_size_ * Vectorizer::SizeOf(instances->nodecount);
    int pdata_size =
        mech->pdata_size_ * Vectorizer::SizeOf(instances->nodecount);

    vector<double> data_new(data_size, 0);
    vector<int> pdata_new(pdata_size, -1);
    vector<int> nodeindices_new(instances->nodecount);
    ml_no_capacitors[m].nodecount = 0;

    // first non-capacitors' instances, then capacitors
    for (int phase = 1; phase <= 2; phase++) {
      for (int n = 0; n < instances->nodecount; n++) {
        int node_id = instances->nodeindices[n];

        // place first the no-caps of valid mechs; then all others
        bool is_capacitor = capacitor_ids.find(node_id) != capacitor_ids.end();
        int instance_phase = !is_capacitor && mech_valid_in_phase_1 ? 1 : 2;

        if (instance_phase != phase) continue;

        assert(n_new < instances->nodecount);
        for (int i = 0; i < mech->data_size_; i++)  // copy data
        {
#if LAYOUT == 1
          int old_data_offset = mech->data_size_ * n + i;
          int new_data_offset = mech->data_size_ * n_new + i;
#else
          int old_data_offset =
              Vectorizer::SizeOf(instances->nodecount) * i + n;
          int new_data_offset =
              Vectorizer::SizeOf(instances->nodecount) * i + n_new;
#endif
          assert(new_data_offset < data_size);
          assert(total_data_offset + old_data_offset ==
                 (&instances->data[old_data_offset] - branch->nt_->_data));
          assert(data_new.at(new_data_offset) == 0);

          data_new.at(new_data_offset) = instances->data[old_data_offset];

          if (mech->is_ion_)
            ions_data_map[mech->type_][total_data_offset + old_data_offset] =
                total_data_offset + new_data_offset;
        }

        for (int i = 0; i < mech->pdata_size_; i++)  // copy pdata
        {
#if LAYOUT == 1
          int old_pdata_offset = mech->pdata_size_ * n + i;
          int new_pdata_offset = mech->pdata_size_ * n_new + i;
#else
          int old_pdata_offset =
              Vectorizer::SizeOf(instances->nodecount) * i + n;
          int new_pdata_offset =
              Vectorizer::SizeOf(instances->nodecount) * i + n_new;
#endif
          int old_pdata = instances->pdata[old_pdata_offset];
          assert(old_pdata >= 0);

          // if it points to an ion, get new pdata position
          int ptype = memb_func[mech->type_].dparam_semantics[i];
          assert(pdata_new.at(new_pdata_offset) == -1);
          if (ptype > 0 && ptype < 1000)  // ptype is ion id
            pdata_new.at(new_pdata_offset) =
                ions_data_map.at(ptype).at(old_pdata);
          else
            pdata_new.at(new_pdata_offset) = old_pdata;
        }

        if (phase == 1)  // count no-caps
          ml_no_capacitors[m].nodecount++;

        nodeindices_new[n_new++] = node_id;
      }
    }
    assert(n_new == nodeindices_new.size());

    // overwite old values in NrnThread->data
    memcpy(ml_no_capacitors[m].data, data_new.data(),
           sizeof(double) * data_new.size());
    memcpy(ml_no_capacitors[m].pdata, pdata_new.data(),
           sizeof(int) * pdata_new.size());
    memcpy(ml_no_capacitors[m].nodeindices, nodeindices_new.data(),
           sizeof(int) * nodeindices_new.size());
    /* Important: do not change padded counts: they are used to properly
     * compute mechanisms offsets (ppvar[x*STRIDE]), ie iterates from 0
     * to count, with gap between mechs variables given by padded_count*/
    total_data_offset += data_size;

    // create Memb_list for capacitors only mechanisms
    if (ml_capacitors_ptr != nullptr) {
      Memb_list*& ml_capacitors = *ml_capacitors_ptr;

      ml_capacitors[m].nodecount =
          instances[m].nodecount - ml_no_capacitors[m].nodecount;
#if LAYOUT == 1
      // address of end of no_cap data
      // AoS data as |abababab| becomes |ababab|+|abab|
      int capacitors_gap_pdata =
          ml_no_capacitors[m].nodecount * mech->pdata_size_;
      int capacitors_gap_data =
          ml_no_capacitors[m].nodecount * mech->data_size_;
#else
      // offset initial address by few positions
      // SoA data as |aaaaabbbbb| becomes |aaa__bbb__|+|__aa__bb|
      int capacitors_gap_data = ml_no_capacitors[m].nodecount;
      int capacitors_gap_pdata = ml_no_capacitors[m].nodecount;
#endif
      ml_capacitors[m].data = &(ml_no_capacitors[m].data[capacitors_gap_data]);
      ml_capacitors[m].pdata =
          &(ml_no_capacitors[m].pdata[capacitors_gap_pdata]);
    }
  }
  assert(total_data_offset == branch->nt_->_ndata);

  // pass output value if required by used
  if (ml_no_capacitors_ptr != nullptr)
    *ml_no_capacitors_ptr = ml_no_capacitors;
  else
    delete[] ml_no_capacitors;
}
