/*
# =============================================================================
# Copyright (c) 2015 - 2021 Blue Brain Project/EPFL
#
# See top-level LICENSE file for details.
# =============================================================================.
*/
#include "neurox/neurox.h"

#include <cmath>
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
      for (size_t i = 0; i < mechanisms_[m]->data_size_; i++)  // and varable
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
        if (input_params_->branch_parallelism_) {
          int ptype = memb_func[mechanisms_[m]->type_].dparam_semantics[i];
          bool is_pointer = ptype == -1 || (ptype > 0 && ptype < 1000);
          if (is_pointer)  // true for pointer to area in nt->data, or ion
                           // instance data
          {
            // if mech has no instances, area is always -1
            if (ptype == -1 /*area*/ &&
                input::DataLoader::HardCodedMechanismHasNoInstances(
                    mechanisms_[m]->type_)) {
              assert(pdata_old[old_offset] == -1);
              pdata_new[new_offset] = pdata_old[old_offset];
            } else {
              // point to new offset
              pdata_new[new_offset] = data_offsets.at(pdata_old[old_offset]);
              assert(data_new[pdata_new[new_offset]] ==
                     data_old[pdata_old[old_offset]]);
            }
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
    int data_size = mech->data_size_ * Vectorizer::SizeOf(instances->nodecount);

    // neuron: "only point processes with currents are possibilities"
    if (input::DataLoader::HardCodedMechanismHasNoInstances(mech->type_)) {
      ml_no_capacitors[m].nodecount = 0;
      ml_no_capacitors[m]._nodecount_padded = 0;
      ml_no_capacitors[m].nodeindices = nullptr;
      ml_no_capacitors[m].data = nullptr;
      ml_no_capacitors[m].nodeindices = nullptr;
      ml_no_capacitors[m].pdata = nullptr;
      total_data_offset += data_size;
      continue;
    }

    bool mech_valid_in_phase_1 = mech->pnt_map_ && mech->memb_func_.current;

    int n_new = 0;
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
          assert(old_pdata >= -1);

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

void tools::Vectorizer::CreateMechInstancesThreads(Branch* b) {
  // First step: measuse compute times
  NrnThread* nt = b->nt_;
  const int benchmark_iterations = 200;
  hpx_time_t now;

  // memory for runtime of individual instances
  std::vector<double> state_func_runtime(mechanisms_count_);
  std::vector<double> current_func_runtime(mechanisms_count_);
  double total_state_instances_runtime = 0;
  double total_current_instances_runtime = 0;

  for (int m = 0; m < neurox::mechanisms_count_; m++) {
    Mechanism* mech = neurox::mechanisms_[m];
    Memb_list* ml = &b->mechs_instances_[m];
    int type = mech->type_;

    /* benchmark state function */
    if (mech->memb_func_.state) {
      now = hpx_time_now();
      for (int i = 0; i < benchmark_iterations; i++)
        mech->memb_func_.state(nt, ml, type);

      // average time per instance per node
      state_func_runtime[m] =
          hpx_time_elapsed_ns(now) / (double)benchmark_iterations;
      total_state_instances_runtime += state_func_runtime[m];
      state_func_runtime[m] /= (double)ml->nodecount;
    }

    /* benchmark current function (capacitance is an exception) */
    if (type == MechanismTypes::kCapacitance || mech->memb_func_.current) {
      now = hpx_time_now();
      for (int i = 0; i < benchmark_iterations; i++)
        mech->memb_func_.current(nt, ml, type, nullptr, nullptr, nullptr);
      current_func_runtime[m] =
          hpx_time_elapsed_ns(now) / (double)benchmark_iterations;
      total_current_instances_runtime += current_func_runtime[m];
      current_func_runtime[m] /= (double)ml->nodecount;
    }
  }

#ifndef NDEBUG
  printf("neuron gid %d: total runtime: state %.3f ns, current %.3f ns\n",
         nt->id, total_state_instances_runtime,
         total_current_instances_runtime);
  for (int m = 0; m < neurox::mechanisms_count_; m++)
    printf("neuron gid %d: mech type %d: state %.8f ns, current %.8f ns\n",
           nt->id, mechanisms_[m]->type_, state_func_runtime[m],
           current_func_runtime[m]);
#endif

  /* compute mech-instances per thread for state and current funtions */
  assert(total_state_instances_runtime > 0 &&
         total_current_instances_runtime > 0);
  const int parallel_simd_instructions_count =
      input_params_->processor_cache_line_size_l1_ / sizeof(double);

  // build memb_lists for parallel processing of mechanisms
  b->mechs_instances_parallel_ =
      new Mechanism::MembListThreadArgs[mechanisms_count_];

  enum funcs { State = 0, Current = 1, Count = 2 };
  for (int f = 0; f < funcs::Count; f++) {
    double max_workload = LoadBalancing::GetWorkloadPerMechInstancesBlock(
        f == funcs::State ? total_state_instances_runtime
                          : total_current_instances_runtime);
    for (int m = 0; m < neurox::mechanisms_count_; m++) {
      Mechanism* mech = neurox::mechanisms_[m];
      Memb_list* ml = &b->mechs_instances_[m];

      Mechanism::MembListThreadArgs* threads_args =
          &b->mechs_instances_parallel_[m];

      /* compute the number of instances per thread thread */
      double instance_runtime =
          f == funcs::State ? state_func_runtime[m] : current_func_runtime[m];

      int cluster_size = ml->nodecount, cluster_count = 1;
      /* if instance_runtime==0, then mod-function is not defined */
      if (instance_runtime > 0 && ml->nodecount > 0) {
        // cluster size is the number of instances that fits the max workload
        //(that fully utilize the processor cache line size)
        cluster_size = std::ceil(max_workload / instance_runtime);
        assert(cluster_size > 0);
        if (cluster_size % parallel_simd_instructions_count != 0)
          cluster_size += parallel_simd_instructions_count -
                          cluster_size % parallel_simd_instructions_count;
        assert(cluster_size % parallel_simd_instructions_count == 0);

        if (cluster_size > ml->nodecount) cluster_size = ml->nodecount;

        cluster_count = std::ceil((double)ml->nodecount / (double)cluster_size);

        // cluster count is the number of parallel threads to be spawned
        assert(cluster_count > 0);
      }
      //#ifndef NDEBUG
      assert(ml->nodecount <= cluster_size * cluster_count);
      assert(ml->nodecount >= cluster_size * (cluster_count - 1));
      printf(
          "neuron gid %d, mech %d, %s: node_count %d, cluster_size %d, "
          "cluster_count %d\n",
          b->soma_->gid_, mech->type_, f == funcs::State ? "state" : "current",
          ml->nodecount, cluster_size, cluster_count);
      //#endif

      /* if cluster_count==1, does not need threaded execution */
      if (instance_runtime == 0 || cluster_count == 1 || ml->nodecount == 0) {
        switch (f) {
          case funcs::State:
            threads_args->ml_state_count = 0;
            threads_args->ml_state = nullptr;
            continue;
            break;
          case funcs::Current:
            threads_args->ml_current_count = 0;
            threads_args->ml_current = nullptr;
            continue;
            break;
          default:
            assert(0);
            break;
        }
      }

      /* if different instances of same mechanism type can be in the same
       * compartment, then instances-parallelism requires mechs-graph parallism
         so that updates to nt->V and nt->RHS are protected by shadow vector */
      if (!input_params_->graph_mechs_parallelism_) {
        Memb_list* ml = &b->mechs_instances_[m];
        int thread_id = 0;

        for (int n = 0; n < ml->nodecount; n += cluster_size, thread_id++) {
          /* map of compartment index to the thread id it's allocated to */
          std::map<int, int> index_to_thread;
          for (int i = 0; i < std::min(ml->nodecount - n, cluster_size); i++) {
            int index = ml->nodeindices[n + i];

            /* if that compartment index can be processed by other thread*/
            if (index_to_thread.find(index) != index_to_thread.end() &&
                index_to_thread.at(index) != thread_id) {
              throw std::runtime_error(
                  "Several instances of same mechanism type in the same "
                  "compartment. Re-run with mechs-graph parallelism.\n");
            } else
              index_to_thread[index] = thread_id;
          }
        }
      }

      threads_args->nt = b->nt_;
      threads_args->memb_func = &mech->memb_func_;
      threads_args->mech_type = mech->type_;

      // for shadow vectors processing of current
      threads_args->requires_shadow_vectors =
          /* graph-parallelism */
          b->mechs_graph_
          /* not CaDynamics_E2 (no updates in cur function) */
          && mech->type_ != MechanismTypes::kCaDynamics_E2
          /* not ion (updates in nrn_cur_ion function) */
          && !mech->is_ion_
          /* not capacitance-only current function */
          &&
          !(mech->type_ == MechanismTypes::kCapacitance && f == funcs::Current);

      /* if graph-parallelism, pass accumulation functions and their argument*/
      if (threads_args->requires_shadow_vectors) {
        threads_args->acc_args = b->mechs_graph_;
        threads_args->acc_rhs_d = Branch::MechanismsGraph::AccumulateRHSandD;

        /* every mech with dependencies will update di_dv of parent mech*/
        threads_args->acc_di_dv =
            mech->dependencies_count_ == 0
                ? nullptr
                : Branch::MechanismsGraph::AccumulateIandDIDV;
      } else {
        threads_args->acc_args = nullptr;
        threads_args->acc_di_dv = nullptr;
        threads_args->acc_rhs_d = nullptr;
      }

      /* parallel subsets of Memb_list */
      int& ml_thread_count = f == funcs::State ? threads_args->ml_state_count
                                               : threads_args->ml_current_count;
      Memb_list*& ml_threads =
          f == funcs::State ? threads_args->ml_state : threads_args->ml_current;

      ml_thread_count = cluster_count;
      ml_threads = new Memb_list[ml_thread_count];

      int thread_id = 0;
      for (int n = 0; n < ml->nodecount; n += cluster_size, thread_id++) {
        Memb_list* thread = &ml_threads[thread_id];
        memcpy(thread, ml, sizeof(Memb_list));
        thread->nodeindices = &ml->nodeindices[n];
        thread->nodecount = std::min(cluster_size, ml->nodecount - n);
        if (input_params_->graph_mechs_parallelism_) {
          thread->_shadow_d = &ml->_shadow_d[n];
          thread->_shadow_rhs = &ml->_shadow_rhs[n];
          thread->_shadow_didv = &ml->_shadow_didv[n];
          thread->_shadow_didv_offsets = &ml->_shadow_didv_offsets[n];
          thread->_shadow_i = &ml->_shadow_i[n];
          thread->_shadow_i_offsets = &ml->_shadow_i_offsets[n];
        }
#if LAYOUT == 1
        int data_offset = n * mech->data_size_;
        int pdata_offset = n * mech->pdata_size_;
#else
        int data_offset = n;
        int pdata_offset = n;
#endif
        thread->data = &ml->data[data_offset];
        thread->pdata = &ml->pdata[pdata_offset];
      }
    }
  }
}
