#include "neurox/neurox.h"

#include <algorithm>
#include <cstring>
#include <numeric>
#include <set>

using namespace neurox;
using namespace neurox::tools;
using namespace neurox::synchronizers;
using namespace neurox::interpolators;

void *Branch::operator new(size_t bytes, void *addr) { return addr; }

void Branch::operator delete(void *worker) {}

Branch::Branch(offset_t n, int nrn_thread_id, int threshold_v_offset,
               floble_t *data, size_t data_count, offset_t *pdata,
               size_t pdata_count, offset_t *instances_count,
               size_t recv_mechs_count, offset_t *nodes_indices,
               size_t nodes_indices_count, hpx_t top_branch_addr,
               hpx_t *branches, size_t branches_count, offset_t *p,
               size_t p_count, floble_t *vecplay_t, size_t vecplay_t_Count,
               floble_t *vecplay_y, size_t vecplay_y_count,
               PointProcInfo *vecplay_ppi, size_t vecplay_ppi_count,
               NetconX *netcons, size_t netcons_count,
               neuron_id_t *netcons_pre_ids, size_t netcons_pre_ids_count,
               floble_t *weights, size_t weights_count,
               unsigned char *vdata_serialized, size_t vdata_serialized_count)
    : soma_(nullptr),
      nt_(nullptr),
      mechs_instances_(nullptr),
      mechs_instances_parallel_(nullptr),
      thvar_ptr_(nullptr),
      mechs_graph_(nullptr),
      branch_tree_(nullptr),
      events_queue_mutex_(HPX_NULL),
      interpolator_(nullptr) {
  this->nt_ = (NrnThread *)malloc(sizeof(NrnThread));
  NrnThread *nt = this->nt_;

  // all non usable values
  nt->_ml_list = NULL;
  nt->tml = NULL;
  nt->pntprocs = NULL;
  nt->presyns = NULL;
  nt->presyns_helper = NULL;
  nt->pnt2presyn_ix = NULL;
  nt->netcons = NULL;
  nt->n_pntproc = -1;
  nt->n_presyn = -1;
  nt->n_input_presyn = -1;
  nt->n_netcon = -1;
  nt->ncell = -1;
  nt->id = nrn_thread_id;
  nt->_stop_stepping = -1;
  nt->_permute = NULL;
  nt->_sp13mat = NULL;
  nt->_ecell_memb_list = NULL;
  nt->_ctime = -1;
  for (int i = 0; i < BEFORE_AFTER_SIZE; i++) nt->tbl[i] = NULL;
  nt->compute_gpu = 0;
  nt->_net_send_buffer_size = -1;
  nt->_net_send_buffer_cnt = -1;
  nt->_net_send_buffer = NULL;
  nt->mapping = NULL;
  nt->_idata = NULL;
  nt->_nidata = -1;

  // assignemnts start here
  nt->_dt = input_params_->dt_;
  nt->_t = input_params_->tstart_;
  nt->cj = (input_params_->second_order_ ? 2.0 : 1.0) / nt->_dt;
  nt->end = n;

  nt->_data = data_count == 0 ? nullptr : Vectorizer::New<floble_t>(data_count);
  memcpy(nt->_data, data, data_count * sizeof(floble_t));
  nt->_ndata = data_count;

  nt->weights = weights_count == 0 ? nullptr : new floble_t[weights_count];
  memcpy(nt->weights, weights, sizeof(floble_t) * weights_count);
  nt->n_weight = weights_count;

  nt->_actual_rhs = nt->_data + n * 0;
  nt->_actual_d = nt->_data + n * 1;
  nt->_actual_a = nt->_data + n * 2;
  nt->_actual_b = nt->_data + n * 3;
  nt->_actual_v = nt->_data + n * 4;
  nt->_actual_area = nt->_data + n * 5;

  // AP threshold offset
  this->thvar_ptr_ = threshold_v_offset == -1
                         ? NULL
                         : &(this->nt_->_actual_v[threshold_v_offset]);

  // events mutex
  this->events_queue_mutex_ = hpx_lco_sema_new(1);

  // parent index
  if (p_count > 0) {
    assert(p_count == n);
    this->nt_->_v_parent_index = Vectorizer::New<offset_t>(p_count);
    memcpy(this->nt_->_v_parent_index, p, n * sizeof(offset_t));
  } else {
    this->nt_->_v_parent_index = nullptr;
  }

  // reconstruct mechanisms
  assert(recv_mechs_count <= mechanisms_count_);
  offset_t data_offset = 6 * n;
  offset_t pdata_offset = 0;
  offset_t instances_offset = 0;
  this->mechs_instances_ = new Memb_list[mechanisms_count_];
  int max_mech_id = 0;

  vector<void *> vdata_ptrs;
  offset_t vdata_offset = 0;

  for (offset_t m = 0; m < mechanisms_count_; m++) {
    Memb_list &instance = this->mechs_instances_[m];
    Mechanism *mech = mechanisms_[m];
    instance.nodecount = instances_count[m];
    max_mech_id = max(max_mech_id, mech->type_);

    // data, pdata, and nodesIndices arrays
    instance.data = mech->data_size_ == 0 || instance.nodecount == 0
                        ? nullptr
                        : this->nt_->_data + data_offset;
    instance.pdata =
        mech->pdata_size_ == 0 || instance.nodecount == 0
            ? nullptr
            : Vectorizer::New<offset_t>(mech->pdata_size_ * instance.nodecount);
    if (instance.pdata)
      memcpy(instance.pdata, &pdata[pdata_offset],
             sizeof(offset_t) * (mech->pdata_size_ * instance.nodecount));
    instance.nodeindices = instance.nodecount > 0
                               ? Vectorizer::New<offset_t>(instance.nodecount)
                               : nullptr;
    if (instance.nodeindices)
      memcpy(instance.nodeindices, &nodes_indices[instances_offset],
             sizeof(offset_t) * instance.nodecount);

    // init thread data :: nrn_setup.cpp->setup_ThreadData();
    if (mech->memb_func_.thread_size_) {
      instance._thread = new ThreadDatum[mech->memb_func_.thread_size_];
      if (mech->memb_func_.thread_mem_init_)
        mech->memb_func_.thread_mem_init_(instance._thread);
    } else {
      instance._thread = nullptr;
    }

    // vdata: if is point process we need to allocate the vdata (by calling
    // bbcore_reg in mod file)
    // and assign the correct offset in pdata (offset of vdata is in pdata[1])
    for (size_t i = 0; i < instance.nodecount; i++) {
      // point pdata to the correct offset, and allocate vdata
      assert(data_offset <= data_count);
      assert(vdata_offset <= vdata_serialized_count);
      floble_t *instance_data = (floble_t *)&this->nt_->_data[data_offset];
      floble_t *instance_data_2 =
          (floble_t *)&instance.data[i * mech->data_size_];
      offset_t *instance_pdata =
          (offset_t *)&instance.pdata[i * mech->pdata_size_];
      assert(instance_data =
                 instance_data_2);  // Make sure data offsets are good so far

      if (mech->vdata_size_ > 0 || mech->pnt_map_ > 0) {
        assert((mech->type_ == MechanismTypes::kIClamp &&
                mech->vdata_size_ == 1 && mech->pdata_size_ == 2 &&
                mech->pnt_map_ > 0) ||
               (mech->type_ == MechanismTypes::kStochKv &&
                mech->vdata_size_ == 1 && mech->pdata_size_ == 5 &&
                mech->pnt_map_ == 0) ||
               ((mech->type_ == MechanismTypes::kProbAMPANMDA_EMS ||
                 mech->type_ == MechanismTypes::kProbGABAAB_EMS) &&
                mech->vdata_size_ == 2 && mech->pdata_size_ == 3 &&
                mech->pnt_map_ > 0));

        // ProbAMPANMDA_EMS, ProbAMPANMDA_EMS and IClamp:
        // pdata[0]: offset in data (area)
        // pdata[1]: offset for Point_process in vdata[0]
        // pdata[2]: offset for RNG in vdata[1]   (NOT for IClamp)

        // StochKv:
        // pdata[0]: offset in area (ion_ek)
        // pdata[1]: offset in area (ion_ik)
        // pdata[2]: offset in area (ion_dikdv)
        // pdata[3]: offset for RNG in vdata[0]
        // pdata[4]: offset in data (area)

        // copy Point_processes by replacing vdata pointers and pdata offset by
        // the ones referring to a copy
        if (mech->type_ == MechanismTypes::kIClamp ||
            mech->type_ == MechanismTypes::kProbAMPANMDA_EMS ||
            mech->type_ == MechanismTypes::kProbGABAAB_EMS) {
          int pointProcOffsetInPdata = 1;
          Point_process *pp =
              (Point_process *)(void *)&vdata_serialized[vdata_offset];
          assert(pp->_i_instance >= 0 && pp->_tid >= 0 && pp->_type >= 0);
          Point_process *pp_copy = new Point_process;
          memcpy(pp_copy, pp, sizeof(Point_process));
          vdata_offset += sizeof(Point_process);
          instance_pdata[pointProcOffsetInPdata] = vdata_ptrs.size();
          vdata_ptrs.push_back(pp_copy);
        }

        // copy RNG by replacing vdata pointers and pdata offset by the ones
        // referring to a copy
        if (mech->type_ == MechanismTypes::kStochKv ||
            mech->type_ == MechanismTypes::kProbAMPANMDA_EMS ||
            mech->type_ == MechanismTypes::kProbGABAAB_EMS) {
          int rng_offset_in_pdata =
              mech->type_ == MechanismTypes::kStochKv ? 3 : 2;
          nrnran123_State *rng =
              (nrnran123_State *)(void *)&vdata_serialized[vdata_offset];
          nrnran123_State *rngcopy = new nrnran123_State;
          memcpy(rngcopy, rng, sizeof(nrnran123_State));
          vdata_offset += sizeof(nrnran123_State);
          instance_pdata[rng_offset_in_pdata] = vdata_ptrs.size();
          vdata_ptrs.push_back(rngcopy);
        }
      }
      data_offset += mech->data_size_;
      pdata_offset += mech->pdata_size_;
      assert(data_offset < 2 ^ sizeof(offset_t));
      assert(pdata_offset < 2 ^ sizeof(offset_t));
      assert(vdata_offset < 2 ^ sizeof(offset_t));
      instances_offset++;
    }
  }
  assert(data_offset == data_count);
  assert(pdata_offset == pdata_count);
  assert(instances_offset == nodes_indices_count);

  // nt->_vdata pointers
  nt->_nvdata = vdata_ptrs.size();
  nt->_vdata =
      vdata_serialized_count == 0 ? nullptr : new void *[vdata_ptrs.size()];
  memcpy(nt->_vdata, vdata_ptrs.data(), vdata_ptrs.size() * sizeof(void *));
  vdata_ptrs.clear();

  // nt->_ml_list
  nt->_ml_list = new Memb_list *[max_mech_id + 1];
  for (int i = 0; i <= max_mech_id; i++) nt->_ml_list[i] = NULL;

  int ionsCount = 0;
  for (offset_t m = 0; m < mechanisms_count_; m++) {
    Mechanism *mech = mechanisms_[m];
    Memb_list &instances = this->mechs_instances_[m];
    this->nt_->_ml_list[mech->type_] = &instances;
    if (mech->is_ion_) ionsCount++;
  }

  // ttx excluded (no writes to ttx state)
  assert(ionsCount == Mechanism::IonTypes::kSizeWriteableIons + 1);

  // vecplay
  nt->n_vecplay = vecplay_ppi_count;
  nt->_vecplay =
      vecplay_ppi_count == 0 ? nullptr : new void *[vecplay_ppi_count];

  offset_t v_offset = 0;
  for (size_t v = 0; v < nt->n_vecplay; v++) {
    PointProcInfo &ppi = vecplay_ppi[v];
    size_t size = ppi.size;
    int m = mechanisms_map_[ppi.mech_type];
    floble_t *instances_data = this->mechs_instances_[m].data;
    assert(mechanisms_[m]->pnt_map_ > 0);
    floble_t *pd =
        &(instances_data[ppi.mech_instance * mechanisms_[m]->data_size_ +
                         ppi.instance_data_offset]);
    floble_t *yvec = new floble_t[size];
    floble_t *tvec = new floble_t[size];
    for (size_t i = 0; i < size; i++) {
      yvec[i] = vecplay_y[v_offset + i];
      tvec[i] = vecplay_t[v_offset + i];
    }
    nt->_vecplay[v] = new VecplayContinuousX(pd, size, yvec, tvec, NULL);
    v_offset += size;
  }

  // Shadow arrays
  this->nt_->shadow_rhs_cnt = 0;
  this->nt_->_shadow_d = NULL;
  this->nt_->_shadow_rhs = NULL;

  for (int m = 0; m < mechanisms_count_; m++) {
    int shadow_size = 0;
    if (mechanisms_[m]->memb_func_.current &&
        !mechanisms_[m]->is_ion_)  // ions have no updates
      shadow_size = this->mechs_instances_[m].nodecount;

    Memb_list *ml = &mechs_instances_[m];
    ml->_shadow_d =
        shadow_size == 0 ? nullptr : Vectorizer::New<double>(shadow_size);
    ml->_shadow_rhs =
        shadow_size == 0 ? nullptr : Vectorizer::New<double>(shadow_size);

    for (int i = 0; i < shadow_size; i++) {
      ml->_shadow_d[i] = 0;
      ml->_shadow_rhs[i] = 0;
    }

    if (mechanisms_[m]->dependency_ion_index_ >=
        Mechanism::IonTypes::kSizeWriteableIons)
      shadow_size = 0;  //> only mechanisms with parent ions update I and DI/DV

    ml->_shadow_i =
        shadow_size == 0 ? nullptr : Vectorizer::New<double>(shadow_size);
    ml->_shadow_didv =
        shadow_size == 0 ? nullptr : Vectorizer::New<double>(shadow_size);
    ml->_shadow_i_offsets =
        shadow_size == 0 ? nullptr : Vectorizer::New<int>(shadow_size);
    ml->_shadow_didv_offsets =
        shadow_size == 0 ? nullptr : Vectorizer::New<int>(shadow_size);
    for (int i = 0; i < shadow_size; i++) {
      ml->_shadow_i[i] = 0;
      ml->_shadow_didv[i] = 0;
      ml->_shadow_i_offsets[i] = -1;
      ml->_shadow_didv_offsets[i] = -1;
    }
  }

  // reconstructs netcons
  offset_t weights_offset = 0;
  for (offset_t nc = 0; nc < netcons_count; nc++) {
    this->netcons_[netcons_pre_ids[nc]].push_back(
        new NetconX(netcons[nc].mech_type_, netcons[nc].mech_instance_,
                    netcons[nc].delay_, netcons[nc].weight_index_,
                    netcons[nc].weights_count_, netcons[nc].active_));
    assert(weights_offset == netcons[nc].weight_index_);
    weights_offset += netcons[nc].weights_count_;
  }
  assert(weights_count == weights_offset);

  // create data structure that defines branching
  if (input_params_->branch_parallelism_)
    this->branch_tree_ =
        new Branch::BranchTree(top_branch_addr, branches, branches_count);

#if LAYOUT == 0
  // if using vector data structures, convert now
  tools::Vectorizer::ConvertToSOA(this);
#endif

  interpolator_ = Interpolator::New(input_params_->interpolator_);
}

void Branch::ClearMembList(Memb_list *&mechs_instances) {
  for (int m = 0; m < mechanisms_count_; m++) {
    Memb_list &instance = mechs_instances[m];
    if (mechanisms_[m]->memb_func_.thread_cleanup_)
      mechanisms_[m]->memb_func_.thread_cleanup_(instance._thread);

    Vectorizer::Delete(mechs_instances[m].nodeindices);
    Vectorizer::Delete(mechs_instances[m].pdata);
    delete[] mechs_instances[m]._thread;

    Vectorizer::Delete(mechs_instances[m]._shadow_d);
    Vectorizer::Delete(mechs_instances[m]._shadow_didv);
    Vectorizer::Delete(mechs_instances[m]._shadow_didv_offsets);
    Vectorizer::Delete(mechs_instances[m]._shadow_i);
    Vectorizer::Delete(mechs_instances[m]._shadow_rhs);
    Vectorizer::Delete(mechs_instances[m]._shadow_i_offsets);
  }
  delete[] mechs_instances;
  mechs_instances = nullptr;
}

void Branch::ClearNrnThread(NrnThread *&nt) {
  delete[] nt->weights;
  Vectorizer::Delete(nt->_data);
  Vectorizer::Delete(nt->_v_parent_index);
  delete[] nt->_ml_list;
  delete[] nt->_vdata;

  for (int i = 0; i < nt->n_vecplay; i++)
    delete (VecplayContinuousX *)nt->_vecplay[i];
  delete[] nt->_vecplay;

  free(nt);
}

Branch::~Branch() {
  hpx_lco_delete_sync(this->events_queue_mutex_);
  ClearMembList(this->mechs_instances_);
  ClearNrnThread(this->nt_);

  for (auto &nc_pair : this->netcons_)
    for (auto &nc : nc_pair.second) delete nc;

  delete soma_;
  delete branch_tree_;
  delete mechs_graph_;
  delete interpolator_;

  if (mechs_instances_parallel_) {
    for (int m = 0; m < neurox::mechanisms_count_; m++) {
      delete[] mechs_instances_parallel_[m].ml_state;
      delete[] mechs_instances_parallel_[m].ml_current;
    }
    delete[] mechs_instances_parallel_;
  }
}

hpx_action_t Branch::Init = 0;
int Branch::Init_handler(const int nargs, const void *args[],
                         const size_t sizes[]) {
  NEUROX_MEM_PIN(Branch);
  assert(nargs == 17);

  new (local) Branch(
      *(offset_t *)args[0],  // number of compartments
      *(int *)args[1],       // nrnThreadId (nt.id)
      *(int *)args[2],       // offset  AP voltage threshold (-1 if none)
      (floble_t *)args[3],
      sizes[3] / sizeof(floble_t),  // data: RHS, D, A, V, B, area, mechs
      (offset_t *)args[4],
      sizes[4] / sizeof(offset_t),  // pdata
      (offset_t *)args[5],
      sizes[5] / sizeof(offset_t),  // instances count per mechanism
      (offset_t *)args[6], sizes[6] / sizeof(offset_t),    // nodes indices
      *(hpx_t *)args[7],                                   // top branchAddr
      (hpx_t *)args[8], sizes[8] / sizeof(hpx_t),          // branches
      (offset_t *)args[9], sizes[9] / sizeof(offset_t),    // parent index
      (floble_t *)args[10], sizes[10] / sizeof(floble_t),  // vecplay T data
      (floble_t *)args[11], sizes[11] / sizeof(floble_t),  // vecplay Y data
      (PointProcInfo *)args[12],
      sizes[12] / sizeof(PointProcInfo),                 // point processes info
      (NetconX *)args[13], sizes[13] / sizeof(NetconX),  // netcons
      (neuron_id_t *)args[14],
      sizes[14] / sizeof(neuron_id_t),                     // netcons pre-ids
      (floble_t *)args[15], sizes[15] / sizeof(floble_t),  // netcons weights
      (unsigned char *)args[16],
      sizes[16] / sizeof(unsigned char));  // serialized vdata
  NEUROX_MEM_UNPIN
}

hpx_action_t Branch::InitSoma = 0;
int Branch::InitSoma_handler(const int nargs, const void *args[],
                             const size_t[]) {
  NEUROX_MEM_PIN(Branch);
  assert(nargs == 2);
  const neuron_id_t neuron_id = *(const neuron_id_t *)args[0];
  const floble_t ap_threshold = *(const floble_t *)args[1];
  local->soma_ = new Neuron(neuron_id, ap_threshold);
  NEUROX_MEM_UNPIN
}

hpx_action_t Branch::Clear = 0;
int Branch::Clear_handler() {
  NEUROX_MEM_PIN(Branch);
  NEUROX_RECURSIVE_BRANCH_ASYNC_CALL(Branch::Clear);
  delete local;
  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
  NEUROX_MEM_UNPIN
}

void Branch::InitVecPlayContinous() {
  // nrn_play_init
  for (size_t v = 0; v < this->nt_->n_vecplay; v++) {
    VecplayContinuousX *vecplay =
        reinterpret_cast<VecplayContinuousX *>(this->nt_->_vecplay[v]);
    vecplay->PlayInit(this);
  }
}

void Branch::AddEventToQueue(floble_t tt, Event *e) {
  this->events_queue_.push(make_pair(tt, e));
}

void Branch::CallModFunction(const Mechanism::ModFunctions function_id,
                             Memb_list *other_ml) {
  if (function_id < BEFORE_AFTER_SIZE) return;  // N/A

  // only for capacitance mechanism
  if (function_id == Mechanism::ModFunctions::kCurrentCapacitance ||
      function_id == Mechanism::ModFunctions::kJacobCapacitance ||
      function_id == Mechanism::ModFunctions::kMulCapacity ||
      function_id == Mechanism::ModFunctions::kDivCapacity) {
    mechanisms_[mechanisms_map_[CAP]]->CallModFunction(this, function_id,
                                                       other_ml);
  }
  // for all others except capacitance (mechanisms graph)
  else {
    if (this->mechs_graph_ != NULL)  // parallel
    {
      // launch execution on top nodes of the branch
      for (int m = 0; m < neurox::mechanisms_count_; m++) {
        if (mechanisms_[m]->type_ == CAP) continue;  // not capacitance
        if (mechanisms_[m]->dependencies_count_ > 0)
          continue;  // not a top branch
        hpx_lco_set(this->mechs_graph_->mechs_lcos_[m], sizeof(function_id),
                    &function_id, HPX_NULL, HPX_NULL);
      }
      // wait for the completion of the graph by waiting at 'end node' lco
      hpx_lco_wait_reset(this->mechs_graph_->end_lco_);
    } else  // serial
    {
      for (int m = 0; m < mechanisms_count_; m++) {
        if (mechanisms_[m]->type_ == CAP &&
            (function_id == Mechanism::ModFunctions::kCurrent ||
             function_id == Mechanism::ModFunctions::kJacob))
          continue;  // kCurrentCapacitance and kJacobCapacitance above
        mechanisms_[m]->CallModFunction(this, function_id, other_ml);
      }
    }
  }
}

hpx_action_t Branch::AddSpikeEventLocality = 0;
int Branch::AddSpikeEventLocality_handler(const int nargs, const void *args[],
                                          const size_t sizes[]) {
  NEUROX_MEM_PIN(uint64_t);
  const neuron_id_t pre_neuron_id = *(const neuron_id_t *)args[0];
  vector<hpx_t> &branch_addrs =
      neurox::locality::netcons_branches_->at(pre_neuron_id);
  hpx_t spikes_lco = hpx_lco_and_new(branch_addrs.size());
  for (hpx_t &branch_addr : branch_addrs)
    if (nargs == 2)
      hpx_call(branch_addr, Branch::AddSpikeEvent, spikes_lco, args[0],
               sizes[0], args[1], sizes[1]);
    else
      hpx_call(branch_addr, Branch::AddSpikeEvent, spikes_lco, args[0],
               sizes[0], args[1], sizes[1], args[2], sizes[2]);
  hpx_lco_wait(spikes_lco);
  hpx_lco_delete(spikes_lco, HPX_NULL);
  NEUROX_MEM_UNPIN
}

// netcvode.cpp::PreSyn::send() --> NetCvode::bin_event.cpp
hpx_action_t Branch::AddSpikeEvent = 0;
int Branch::AddSpikeEvent_handler(const int nargs, const void *args[],
                                  const size_t[]) {
  NEUROX_MEM_PIN(Branch);
  const neuron_id_t pre_neuron_id = *(const neuron_id_t *)args[0];
  const spike_time_t spike_time = *(const spike_time_t *)args[1];
  spike_time_t dependency_time =
      nargs == 3 ? *(const spike_time_t *)args[2] : -1;

  assert(local->netcons_.find(pre_neuron_id) != local->netcons_.end());
  auto &netcons = local->netcons_.at(pre_neuron_id);
  hpx_lco_sema_p(local->events_queue_mutex_);
  for (auto nc : netcons) {
    floble_t delivery_time = spike_time + nc->delay_;
    local->events_queue_.push(make_pair(delivery_time, (Event *)nc));
  }
  hpx_lco_sema_v_sync(local->events_queue_mutex_);

  synchronizer_->AfterReceiveSpikes(local, target, pre_neuron_id, spike_time,
                                    dependency_time);
  NEUROX_MEM_UNPIN;
}

hpx_action_t Branch::ThreadTableCheck = 0;
int Branch::ThreadTableCheck_handler() {
  NEUROX_MEM_PIN(Branch);
  NEUROX_RECURSIVE_BRANCH_ASYNC_CALL(Branch::ThreadTableCheck);
  local->CallModFunction(Mechanism::ModFunctions::kThreadTableCheck);
  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
  NEUROX_MEM_UNPIN;
}

void Branch::DeliverEvents(floble_t til)  // Coreneuron: til=t+0.5*dt
{
  // delivers events in the previous half-step
  floble_t tsav = this->nt_->_t;  // copying cvodestb.cpp logic
  hpx_lco_sema_p(this->events_queue_mutex_);
  while (!this->events_queue_.empty() &&
         this->events_queue_.top().first <= til) {
    auto event_it = this->events_queue_.top();
    floble_t &tt = event_it.first;
    Event *&e = event_it.second;

    // must have been delivered in the previous half step (not earlier), or it
    // was missed
    // TODO assert(tt >= til-0.5*nt->_dt && tt<=til);

    e->Deliver(tt, this);
    this->events_queue_.pop();
  }
  hpx_lco_sema_v_sync(this->events_queue_mutex_);
  this->nt_->_t = tsav;
}

void Branch::FixedPlayContinuous(double t) {
  for (int v = 0; v < this->nt_->n_vecplay; v++) {
    void *vecplay_void = this->nt_->_vecplay[v];
    VecplayContinuousX *vecplay =
        reinterpret_cast<VecplayContinuousX *>(vecplay_void);
    vecplay->Continuous(t);
  }
}

void Branch::FixedPlayContinuous() { FixedPlayContinuous(this->nt_->_t); }

hpx_action_t Branch::SetSyncStepTrigger = 0;
int Branch::SetSyncStepTrigger_handler(const hpx_t *step_trigger_ptr,
                                       const size_t) {
  NEUROX_MEM_PIN(Branch);
  assert(local->soma_);
  local->soma_->synchronizer_step_trigger_ = *step_trigger_ptr;
  NEUROX_MEM_UNPIN
}

//////////////////// Branch::NeuronTree ///////////////////////

Branch::BranchTree::BranchTree(hpx_t top_branch_addr, hpx_t *branches,
                               size_t branches_count)
    : top_branch_addr_(top_branch_addr),
      branches_(nullptr),
      branches_count_(branches_count),
      with_children_lcos_(nullptr),
      children_a_(nullptr),
      children_b_(nullptr),
      parent_a_(-1),
      parent_b_(-1) {
  if (branches_count > 0) {
    this->branches_ = new hpx_t[branches_count];
    memcpy(this->branches_, branches, branches_count * sizeof(hpx_t));
  }
}

Branch::BranchTree::~BranchTree() {
  delete[] branches_;
  delete[] with_children_lcos_;
  delete[] children_a_;
}

hpx_action_t Branch::InitMechanismsGraph = 0;
int Branch::InitMechanismsGraph_handler() {
  NEUROX_MEM_PIN(Branch);
  NEUROX_RECURSIVE_BRANCH_ASYNC_CALL(Branch::InitMechanismsGraph);
  // create data structure that defines the graph of mechanisms
  if (input_params_->graph_mechs_parallelism_) {
    local->mechs_graph_ = new Branch::MechanismsGraph();
    for (size_t m = 0; m < mechanisms_count_; m++) {
      if (mechanisms_[m]->type_ == CAP) continue;  // exclude capacitance
      hpx_call(target, MechanismsGraph::MechFunction,
               local->mechs_graph_->graph_lco_, &mechanisms_[m]->type_,
               sizeof(int));
    }
  }
  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
  NEUROX_MEM_UNPIN;
};

hpx_action_t Branch::InitMechParallelism = 0;
int Branch::InitMechParallelism_handler() {
  NEUROX_MEM_PIN(Branch);
  NEUROX_RECURSIVE_BRANCH_ASYNC_CALL(Branch::InitMechParallelism);
  // creates data structures that defines mech-instances parallelism
  // (requires data structures to be final, so it's run as last)
  if (input_params_->mech_instances_parallelism_)
    Vectorizer::CreateMechInstancesThreads(local);
  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
  NEUROX_MEM_UNPIN;
};

hpx_action_t Branch::Initialize = 0;
int Branch::Initialize_handler() {
  NEUROX_MEM_PIN(Branch);
  NEUROX_RECURSIVE_BRANCH_ASYNC_CALL(Branch::Initialize);
  HinesSolver::CommunicateConstants(local);
  local->interpolator_->Init(local);  // Finitialize, Cvodes init, etc
  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
  NEUROX_MEM_UNPIN;
}

hpx_action_t Branch::BranchTree::InitLCOs = 0;
int Branch::BranchTree::InitLCOs_handler() {
  NEUROX_MEM_PIN(Branch);
  BranchTree *branch_tree = local->branch_tree_;
  if (branch_tree) {
    offset_t branches_count = branch_tree->branches_count_;
    for (int i = 0; i < BranchTree::kFuturesSize; i++) {
      branch_tree->with_parent_lco_[i] = HPX_NULL;
      if (!local->soma_)  // create the LCOs
      {
        size_t size_buffer = i >= 3 ? 2 : 1;
        branch_tree->with_parent_lco_[i] =
            hpx_lco_future_new(size_buffer * sizeof(floble_t));
      }
    }
    branch_tree->with_children_lcos_ =
        branches_count ? new hpx_t[branches_count][BranchTree::kFuturesSize]
                       : nullptr;

    // send my LCOs to children, and receive theirs
    if (branches_count > 0) {
      hpx_t *futures = branches_count ? new hpx_t[branches_count] : nullptr;
      void **addrs = branches_count ? new void *[branches_count] : nullptr;
      size_t *sizes = branches_count ? new size_t[branches_count] : nullptr;
      for (offset_t c = 0; c < branches_count; c++) {
        futures[c] =
            hpx_lco_future_new(sizeof(hpx_t) * BranchTree::kFuturesSize);
        addrs[c] = &branch_tree->with_children_lcos_[c];
        sizes[c] = sizeof(hpx_t) * BranchTree::kFuturesSize;
        hpx_call(
            branch_tree->branches_[c], Branch::BranchTree::InitLCOs, futures[c],
            branch_tree->with_parent_lco_,
            sizeof(hpx_t) * BranchTree::kFuturesSize);  // pass my LCOs down
      }
      hpx_lco_get_all(branches_count, futures, sizes, addrs, nullptr);
      hpx_lco_delete_all(branches_count, futures, HPX_NULL);

      delete[] futures;
      delete[] addrs;
      delete[] sizes;
    }

    if (!local->soma_)  // send my LCO to parent
      return neurox::wrappers::MemoryUnpin(target,
                                           branch_tree->with_parent_lco_);
  }
  NEUROX_MEM_UNPIN;
}

//////////////////// Branch::MechanismsGraph ///////////////////////

Branch::MechanismsGraph::MechanismsGraph() {
  // initializes mechanisms graphs (capacitance is excluded from graph)
  this->graph_lco_ =
      hpx_lco_and_new(mechanisms_count_ - 1);  // excludes 'capacitance'
  this->mechs_lcos_ = new hpx_t[mechanisms_count_];
  size_t terminal_mechanisms_count = 0;
  for (size_t m = 0; m < mechanisms_count_; m++) {
    this->mechs_lcos_[m] = HPX_NULL;
    if (mechanisms_[m]->type_ == CAP) continue;  // exclude capacitance

    this->mechs_lcos_[m] = hpx_lco_reduce_new(
        max((short)1, mechanisms_[m]->dependencies_count_),
        sizeof(Mechanism::ModFunctions), Branch::MechanismsGraph::Init,
        Branch::MechanismsGraph::Reduce);
    if (mechanisms_[m]->successors_count_ == 0)  // bottom of mechs graph
      terminal_mechanisms_count++;
  }
  this->end_lco_ = hpx_lco_and_new(terminal_mechanisms_count);

  this->rhs_d_mutex_ = hpx_lco_sema_new(1);
  for (int i = 0; i < Mechanism::IonTypes::kSizeWriteableIons; i++)
    this->i_didv_mutex_[i] = hpx_lco_sema_new(1);
}

Branch::MechanismsGraph::~MechanismsGraph() {
  hpx_lco_delete_sync(end_lco_);
  hpx_lco_delete_sync(graph_lco_);

  for (int i = 0; i < mechanisms_count_; i++)
    if (mechs_lcos_[i] != HPX_NULL) hpx_lco_delete_sync(mechs_lcos_[i]);
  delete[] mechs_lcos_;

  hpx_lco_delete_sync(rhs_d_mutex_);
  for (int i = 0; i < Mechanism::IonTypes::kSizeWriteableIons; i++)
    hpx_lco_delete_sync(i_didv_mutex_[i]);
}

hpx_action_t Branch::MechanismsGraph::MechFunction = 0;
int Branch::MechanismsGraph::MechFunction_handler(const int *mech_type_ptr,
                                                  const size_t) {
  NEUROX_MEM_PIN(Branch);
  int type = *mech_type_ptr;
  assert(type != CAP);  // capacitance should be outside mechanisms graph
  assert(local->mechs_graph_);
  assert(local->mechs_graph_->mechs_lcos_[mechanisms_map_[type]] != HPX_NULL);
  Mechanism *mech = GetMechanismFromType(type);

  Mechanism::ModFunctions function_id;
  while (local->mechs_graph_->graph_lco_ != HPX_NULL) {
    // wait until all dependencies have completed, and get the argument
    //(function id) from the hpx_lco_set
    hpx_lco_get_reset(local->mechs_graph_->mechs_lcos_[mechanisms_map_[type]],
                      sizeof(Mechanism::ModFunctions), &function_id);
    assert(function_id != Mechanism::ModFunctions::kJacobCapacitance);
    assert(function_id != Mechanism::ModFunctions::kCurrentCapacitance);
    assert(function_id != Mechanism::ModFunctions::kDivCapacity);
    mech->CallModFunction(local, function_id);

    if (mech->successors_count_ == 0)  // bottom mechanism
      hpx_lco_set(local->mechs_graph_->end_lco_, 0, NULL, HPX_NULL, HPX_NULL);
    else
      for (int c = 0; c < mech->successors_count_; c++)
        hpx_lco_set(local->mechs_graph_
                        ->mechs_lcos_[mechanisms_map_[mech->successors_[c]]],
                    sizeof(function_id), &function_id, HPX_NULL, HPX_NULL);
  }
  NEUROX_MEM_UNPIN;
}

void Branch::MechanismsGraph::AccumulateRHSandD(NrnThread *nt, Memb_list *ml,
                                                int, void *args) {
  MechanismsGraph *mg = (MechanismsGraph *)args;
  hpx_lco_sema_p(mg->rhs_d_mutex_);
  for (int n = 0; n < ml->nodecount; n++) {
    int &i = ml->nodeindices[n];
    nt->_actual_rhs[i] += ml->_shadow_rhs[n];
    nt->_actual_d[i] += ml->_shadow_d[n];
  }
  hpx_lco_sema_v_sync(mg->rhs_d_mutex_);
}

void Branch::MechanismsGraph::AccumulateIandDIDV(NrnThread *nt, Memb_list *ml,
                                                 int type, void *args) {
  MechanismsGraph *mg = (MechanismsGraph *)args;
  Mechanism *mech = GetMechanismFromType(type);
  assert(mech->dependency_ion_index_ < Mechanism::IonTypes::kSizeWriteableIons);
  hpx_lco_sema_p(mg->i_didv_mutex_[mech->dependency_ion_index_]);
  for (int n = 0; n < ml->nodecount; n++) {
    int &i_offset = ml->_shadow_i_offsets[n];
    int &didv_offset = ml->_shadow_didv_offsets[n];
    assert(i_offset >= 0 && didv_offset >= 0);
    nt->_data[i_offset] += ml->_shadow_i[n];
    nt->_data[didv_offset] += ml->_shadow_didv[n];
  }
  hpx_lco_sema_v_sync(mg->i_didv_mutex_[mech->dependency_ion_index_]);
}

hpx_action_t Branch::MechanismsGraph::Init = 0;
void Branch::MechanismsGraph::Init_handler(Mechanism::ModFunctions *,
                                           const size_t) {}

hpx_action_t Branch::MechanismsGraph::Reduce = 0;
void Branch::MechanismsGraph::Reduce_handler(Mechanism::ModFunctions *lhs,
                                             const Mechanism::ModFunctions *rhs,
                                             const size_t) {
  *lhs = *rhs;
}

void Branch::RegisterHpxActions() {
  wrappers::RegisterZeroVarAction(Branch::Clear, Branch::Clear_handler);
  wrappers::RegisterZeroVarAction(Branch::ThreadTableCheck,
                                  Branch::ThreadTableCheck_handler);
  wrappers::RegisterZeroVarAction(BranchTree::InitLCOs,
                                  BranchTree::InitLCOs_handler);
  wrappers::RegisterZeroVarAction(Branch::InitMechanismsGraph,
                                  Branch::InitMechanismsGraph_handler);
  wrappers::RegisterZeroVarAction(Branch::InitMechParallelism,
                                  Branch::InitMechParallelism_handler);
  wrappers::RegisterZeroVarAction(Branch::Initialize,
                                  Branch::Initialize_handler);
  wrappers::RegisterSingleVarAction<int>(MechanismsGraph::MechFunction,
                                         MechanismsGraph::MechFunction_handler);
  wrappers::RegisterSingleVarAction<hpx_t>(Branch::SetSyncStepTrigger,
                                           Branch::SetSyncStepTrigger_handler);
  wrappers::RegisterMultipleVarAction(Branch::Init, Branch::Init_handler);
  wrappers::RegisterMultipleVarAction(Branch::InitSoma,
                                      Branch::InitSoma_handler);
  wrappers::RegisterMultipleVarAction(Branch::AddSpikeEvent,
                                      Branch::AddSpikeEvent_handler);
  wrappers::RegisterMultipleVarAction(Branch::AddSpikeEventLocality,
                                      Branch::AddSpikeEventLocality_handler);
  wrappers::RegisterAllReduceInitAction<Mechanism::ModFunctions>(
      Branch::MechanismsGraph::Init, Branch::MechanismsGraph::Init_handler);
  wrappers::RegisterAllReduceReduceAction<Mechanism::ModFunctions>(
      Branch::MechanismsGraph::Reduce, Branch::MechanismsGraph::Reduce_handler);
}
