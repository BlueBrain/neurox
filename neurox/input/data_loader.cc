#include "neurox/neurox.h"

#include <stdio.h>
#include <algorithm>
#include <list>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <tuple>

#include "coreneuron/nrnconf.h"
#include "coreneuron/nrniv/netcon.h"
#include "coreneuron/nrniv/netcvode.h"
#include "coreneuron/nrniv/nrn_assert.h"
#include "coreneuron/nrniv/nrn_setup.h"
#include "coreneuron/nrniv/nrniv_decl.h"
#include "coreneuron/nrniv/vrecitem.h"
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrnoc/nrnoc_decl.h"  //nrn_is_ion()
#include "coreneuron/utils/memory_utils.h"
#include "coreneuron/utils/randoms/nrnran123.h"  //RNG data structures

using namespace std;
using namespace neurox::input;
using namespace neurox::synchronizers;
using namespace neurox::tools;

FILE *DataLoader::file_netcons_ = nullptr;
hpx_t DataLoader::locality_mutex_ = HPX_NULL;
std::vector<hpx_t> *DataLoader::my_neurons_addr_ = nullptr;
std::vector<int> *DataLoader::my_neurons_gids_ = nullptr;
std::vector<int> *DataLoader::all_neurons_gids_ = nullptr;
tools::LoadBalancing *DataLoader::load_balancing_ = nullptr;

neuron_id_t DataLoader::GetNeuronIdFromNrnThreadId(int nrn_id) {
  return (neuron_id_t)nrn_threads[nrn_id].presyns[0].gid_;
}

void DataLoader::PrintSubClustersToFile(FILE *file_compartments,
                                        Compartment *top_compartment) {
  if (input_params_->output_compartments_dot_) {
    assert(top_compartment != NULL);
    fprintf(file_compartments, "subgraph cluster_%d { ", top_compartment->id_);
    if (top_compartment->id_ == 0)
      fprintf(file_compartments, "label=\"SOMA\"; ");
    else if (top_compartment->id_ == 2)
      fprintf(file_compartments, "label=\"AIS\"; ");
    Compartment *comp = NULL;
    for (comp = top_compartment; comp->branches_.size() == 1;
         comp = comp->branches_.front())
      fprintf(file_compartments, "%d; ", comp->id_);
    fprintf(file_compartments, "%d };\n", comp->id_);
    for (int c = 0; c < comp->branches_.size(); c++)
      PrintSubClustersToFile(file_compartments, comp->branches_[c]);
  }
}

PointProcInfo DataLoader::GetPointProcInfoFromDataPointer(NrnThread *nt,
                                                          double *pd,
                                                          size_t size) {
  PointProcInfo ppi;
  ppi.node_id = -1;
  ppi.size = size;
  bool found = false;
  for (NrnThreadMembList *tml = nt->tml; !found && tml != NULL;
       tml = tml->next)  // For every mechanism
  {
    int type = tml->index;
    Mechanism *mech = GetMechanismFromType(type);
    Memb_list *ml = tml->ml;
    for (int n = 0; n < ml->nodecount; n++)
      for (int i = 0; i < mech->data_size_; i++) {
#if LAYOUT == 1
        int data_offset = mech->data_size_ * n + i;
#else
        int data_offset = Vectorizer::SizeOf(ml->nodecount) * i + n;
#endif
        if (&ml->data[data_offset] != pd) continue;  // if not this variable

        ppi.node_id = ml->nodeindices[n];
        ppi.mech_type = type;
        ppi.mech_instance = (offset_t)n;
        ppi.instance_data_offset = i;
        found = true;
        break;
      }
  }
  assert(found);
  return ppi;
}

int DataLoader::CreateNeuron(int neuron_idx, void *) {
  NrnThread *nt = &nrn_threads[neuron_idx];
  neuron_id_t neuron_id = GetNeuronIdFromNrnThreadId(nt->id);
  int N = nt->end;

  // if data is permuted, this method fails.
  assert(!use_interleave_permute && !use_solve_interleave &&
         nt->_permute == NULL);

  // map of padded to non-padded offsets of data
  size_t data_size_padded = 6 * Vectorizer::SizeOf(N);
  for (NrnThreadMembList *tml = nt->tml; tml != NULL; tml = tml->next)
    data_size_padded += Vectorizer::SizeOf(tml->ml->nodecount) *
                        mechanisms_[mechanisms_map_[tml->index]]->data_size_;

  // map of post- to pre-padding values of pdata
  // only used for branched neurons, otherwise pointers for padded and
  // non-padded layouts are the same
  std::vector<int> data_offsets(
      input_params_->branch_parallelism_complexity_ > 0 ? data_size_padded : 0,
      -99999);

  if (input_params_->branch_parallelism_complexity_ > 0)
    for (int n = 0; n < N; n++)
      for (int i = 0; i < 6; i++) {
        int offset_padded = Vectorizer::SizeOf(N) * i + n;
        int offset_non_padded = N * i + n;
        data_offsets[offset_padded] = offset_non_padded;
      }

  //======= 1 - reconstructs matrix (solver) values =======
  deque<Compartment *> compartments;
  for (int n = 0; n < nt->end; n++)
    compartments.push_back(new Compartment(
        (offset_t)n, (floble_t)nt->_actual_a[n], (floble_t)nt->_actual_b[n],
        (floble_t)nt->_actual_d[n], (floble_t)nt->_actual_v[n],
        (floble_t)nt->_actual_rhs[n], (floble_t)nt->_actual_area[n],
        (offset_t)nt->_v_parent_index[n]));

  // reconstructs parents tree
  for (int n = 1; n < N; n++)  // exclude top (no parent)
  {
    Compartment *parent_compartment = compartments.at(nt->_v_parent_index[n]);
    parent_compartment->AddChild(compartments.at(n));
  }

  if (input_params_->output_compartments_dot_) {
    set<int> axon_initial_segment_compartments;
    FILE *file_compartments = fopen(
        string("compartments_" + to_string(neuron_id) + ".dot").c_str(), "wt");
    fprintf(file_compartments, "graph G_%d\n{ bgcolor=%s; \n", neuron_id,
            "transparent");
    fprintf(file_compartments,
            "graph [fontname=helvetica, style=filled, color=blue, "
            "fillcolor=floralwhite];\n");
    fprintf(file_compartments,
            "node [fontname=helvetica, shape=cylinder, color=gray, "
            "style=filled, fillcolor=white];\n");
    fprintf(file_compartments, "edge [fontname=helvetica, color=gray];\n");
    PrintSubClustersToFile(file_compartments,
                           compartments.at(0));  // add subclusters
    for (Compartment *comp : compartments)       // draw edges
      for (int c = 0; c < comp->branches_.size(); c++) {
        bool isSoma = comp->id_ == 1;  // bottom of soma
        Compartment *child = comp->branches_.at(c);
        if ((isSoma && c == 0)  // connection from some to AIS
            ||
            axon_initial_segment_compartments.find(comp->id_) !=
                axon_initial_segment_compartments
                    .end())  // connection to any AIS compartment
        {
          // reverse connection, so that it plots AIS on top of soma in dot file
          fprintf(file_compartments, "%d -- %d%s;\n", child->id_, comp->id_,
                  comp->branches_.size() == 1 ? "" : " [color=black]");
          axon_initial_segment_compartments.insert(child->id_);
        } else
          fprintf(file_compartments, "%d -- %d%s;\n", comp->id_, child->id_,
                  comp->branches_.size() == 1 ? "" : " [color=black]");
      }
    fprintf(file_compartments, "}\n");
    fclose(file_compartments);
  }

  //======= 2 - reconstructs mechanisms instances ========
  unsigned vdata_total_offset = 0;
  unsigned data_total_offset = N * 6;                             // no padding
  unsigned data_total_padded_offset = Vectorizer::SizeOf(N) * 6;  // with
                                                                  // padding
  unsigned point_proc_total_offset = 0;

  // information about offsets in data and node ifs of all instances of all ions
  vector<DataLoader::IonInstancesInfo> ions_instances_info(
      Mechanism::IonTypes::kSizeAllIons);
  for (NrnThreadMembList *tml = nt->tml; tml != NULL;
       tml = tml->next)  // For every mechanism
  {
    int type = tml->index;
    Memb_list *ml = tml->ml;  // Mechanisms application to each compartment
    Mechanism *mech = GetMechanismFromType(type);
    assert(mech->type_ == type);
    int ion_offset = mech->GetIonIndex();

    // for every mech instance (or compartment this mech is applied to)
    for (int n = 0; n < ml->nodecount; n++) {
      if (mech->is_ion_) {
        if (n == 0) {
          ions_instances_info[ion_offset].mech_type = type;
          ions_instances_info[ion_offset].data_start = data_total_offset;
          ions_instances_info[ion_offset].data_end =
              data_total_offset + ml->nodecount * mech->data_size_;
        }
        ions_instances_info[ion_offset].node_ids.push_back(ml->nodeindices[n]);
      }

      std::vector<double> data;
      for (int i = 0; i < mech->data_size_; i++) {
        int offset_non_padded = mech->data_size_ * n + i;
#if LAYOUT == 1
        data.push_back(ml->data[offset_non_padded]);
        assert(ml->data[offset_non_padded] ==
               nt->_data[data_total_offset + offset_non_padded]);
#else
        int offset_padded = Vectorizer::SizeOf(ml->nodecount) * i + n;
        data.push_back(ml->data[offset_padded]);
        assert(ml->data[offset_padded] ==
               nt->_data[data_total_padded_offset + offset_padded]);
        if (input_params_->branch_parallelism_complexity_ > 0)
          data_offsets[data_total_padded_offset + offset_padded] =
              data_total_offset + offset_non_padded;
#endif
      }
      data.shrink_to_fit();

      std::vector<int> pdata;
      for (int i = 0; i < mech->pdata_size_; i++) {
#if LAYOUT == 1
        int pdata_offset_non_padded = mech->pdata_size_ * n + i;
        pdata.push_back(ml->pdata[pdata_offset_non_padded]);
#else
        int pdata_offset_padded = Vectorizer::SizeOf(ml->nodecount) * i + n;
        int pd = ml->pdata[pdata_offset_padded];
        int ptype = memb_func[mech->type_].dparam_semantics[i];

        // remove extra space added by padding (for pointer to area or ion mech
        // instance)
        if (input_params_->branch_parallelism_complexity_ > 0 &&
            (ptype == -1 || (ptype > 0 && ptype < 1000))) {
          assert(data_offsets.at(pd) != -99999);
          pdata.push_back(data_offsets.at(pd));  // offset to non-padded SoA
                                                 // value
        } else
          pdata.push_back(pd);
#endif
      }
      pdata.shrink_to_fit();

      void **vdata = &nt->_vdata[vdata_total_offset];
      Compartment *compartment = compartments.at(ml->nodeindices[n]);
      assert(compartment->id_ == ml->nodeindices[n]);

      if (mech->pnt_map_ > 0 || mech->vdata_size_ > 0)  // vdata
      {
        assert((type == MechanismTypes::kIClamp && mech->vdata_size_ == 1 &&
                mech->pdata_size_ == 2 && mech->pnt_map_ > 0) ||
               (type == MechanismTypes::kStochKv && mech->vdata_size_ == 1 &&
                mech->pdata_size_ == 5 && mech->pnt_map_ == 0) ||
               ((type == MechanismTypes::kProbAMPANMDA_EMS ||
                 type == MechanismTypes::kProbGABAAB_EMS) &&
                mech->vdata_size_ == 2 && mech->pdata_size_ == 3 &&
                mech->pnt_map_ > 0));

        // ProbAMPANMDA_EMS, ProbAMPANMDA_EMS and IClamp:
        // pdata[0]: offset in data (area)
        // pdata[1]: offset for Point_process in vdata[0]
        // pdata[2]: offset for RNG in vdata[1]   (NOT for IClamp,  =pdata[1]+1)

        // StochKv:
        // pdata[0]: offset in area (ion_ek)
        // pdata[1]: offset in area (ion_ik)
        // pdata[2]: offset in area (ion_dikdv)
        // pdata[3]: offset for RNG in vdata[0]
        // pdata[4]: offset in data (area)

        if (type == MechanismTypes::kIClamp ||
            type == MechanismTypes::kProbAMPANMDA_EMS ||
            type == MechanismTypes::kProbGABAAB_EMS) {
          const int point_proc_offset_in_pdata = 1;
          Point_process *ppn = &nt->pntprocs[point_proc_total_offset++];
          assert(nt->_vdata[pdata[point_proc_offset_in_pdata]] == ppn);
          Point_process *pp = (Point_process *)vdata[0];
          assert(pp != NULL);
          compartment->AddSerializedVData((unsigned char *)(void *)pp,
                                          sizeof(Point_process));
        }

        if (type == MechanismTypes::kStochKv ||
            type == MechanismTypes::kProbAMPANMDA_EMS ||
            type == MechanismTypes::kProbGABAAB_EMS) {
          int rng_offset_in_pdata =
              mech->type_ == MechanismTypes::kStochKv ? 3 : 2;
          int rng_offset_in_vdata =
              mech->type_ == MechanismTypes::kStochKv ? 0 : 1;
          assert(nt->_vdata[pdata[rng_offset_in_pdata]] ==
                 vdata[rng_offset_in_vdata]);
          nrnran123_State *RNG = (nrnran123_State *)vdata[rng_offset_in_vdata];

          // TODO: manual hack: StochKv's current state has NULL pointers, why?
          if (RNG == NULL && type == MechanismTypes::kStochKv)
            compartment->AddSerializedVData(
                (unsigned char *)new nrnran123_State, sizeof(nrnran123_State));
          else {
            assert(RNG != NULL);
            compartment->AddSerializedVData((unsigned char *)(void *)RNG,
                                            sizeof(nrnran123_State));
          }
        }
      }
      compartment->AddMechanismInstance(type, n, data.data(), mech->data_size_,
                                        pdata.data(), mech->pdata_size_);

      vdata_total_offset += (unsigned)mech->vdata_size_;
    }
    data_total_offset += mech->data_size_ * ml->nodecount;
    data_total_padded_offset +=
        mech->data_size_ * Vectorizer::SizeOf(ml->nodecount);
  }
  for (Compartment *comp : compartments) comp->ShrinkToFit();

  data_offsets.clear();

  //======= 3 - reconstruct NetCons =====================
  map<neuron_id_t, vector<NetconX *>>
      netcons;  // netcons per pre-synaptic neuron id)
  for (int n = 0; n < nt->n_netcon; ++n) {
    NetCon *nc = nt->netcons + n;
    assert(netcon_srcgid.size() >
           0);  // if size==0, then setup_cleanup() in nrn_setup.cpp was called
    assert(nt->id < netcon_srcgid.size());
    assert(netcon_srcgid.at(nt->id) != NULL);
    int srcgid = netcon_srcgid[nt->id][n];
    assert(srcgid >= 0);  // gids can be negative! this is a reminder that i
                          // should double-check when they happen
    int mech_type = nc->target_->_type;
    int weights_count = pnt_receive_size[mech_type];
    size_t weights_offset = nc->u.weight_index_;
    assert(weights_offset < nt->n_weight);

    int node_id = nrn_threads[nc->target_->_tid]
                      ._ml_list[mech_type]
                      ->nodeindices[nc->target_->_i_instance];
    Compartment *comp = compartments.at(node_id);
    NetconX *nx = new NetconX(mech_type, (offset_t)nc->target_->_i_instance,
                              (floble_t)nc->delay_, weights_offset,
                              weights_count, nc->active_);
    double *weights = &nt->weights[weights_offset];
    comp->AddNetcon(srcgid, nx, weights);
    netcons[srcgid].push_back(nx);
  }

  if (input_params_->output_netcons_dot) {
    int netcons_from_others = 0;
    for (auto nc : netcons) {
      int srcGid = nc.first;
      if (std::find(all_neurons_gids_->begin(), all_neurons_gids_->end(),
                    srcGid) == all_neurons_gids_->end())
        netcons_from_others++;
      else {
        floble_t min_delay = 99999;
        for (auto ncv : nc.second)  // get minimum delay between neurons
          min_delay = std::min(min_delay, ncv->delay_);
        fprintf(file_netcons_, "%d -> %d [label=\"%d (%.2fms)\"];\n", srcGid,
                neuron_id, nc.second.size(), min_delay);
      }
    }
#if NEUROX_INPUT_DATALOADER_OUTPUT_EXTERNAL_NETCONS == true
    if (netcons_from_others > 0)
      fprintf(file_netcons_,
              "%s -> %d [label=\"%d\" fontcolor=gray color=gray arrowhead=vee "
              "fontsize=12];\n",
              "external", neuron_id, netcons_from_others);
#endif
  }

  for (auto nc : netcons) {
    for (auto ncv : nc.second) delete ncv;
    nc.second.clear();
  }
  netcons.clear();
  for (Compartment *comp : compartments) comp->ShrinkToFit();

  //======= 4 - reconstruct VecPlayContinuous events =======
  for (int v = 0; v < nt->n_vecplay; v++) {
    VecPlayContinuous *vec = (VecPlayContinuous *)nt->_vecplay[v];
    // discover node, mechanism and data offset id that *pd points to
    PointProcInfo ppi =
        GetPointProcInfoFromDataPointer(nt, vec->pd_, vec->y_->size());
    compartments.at(ppi.node_id)
        ->AddVecplay(vec->t_->data(), vec->y_->data(), ppi);
  }

  //======= 5 - recursively create branches tree ===========
  floble_t ap_threshold = (floble_t)nrn_threads[nt->id].presyns[0].threshold_;
  int thvar_index = nrn_threads[nt->id].presyns[0].thvar_index_;

  for (Compartment *comp : compartments) comp->ShrinkToFit();

  CreateBranch(nt->id, HPX_NULL, compartments, compartments.at(0),
               ions_instances_info, -1, thvar_index, ap_threshold);

  for (auto c : compartments) delete c;
  for (auto nc : netcons)
    for (auto nvc : nc.second) delete nvc;
  return 0;
}

void DataLoader::CleanCoreneuronData(const bool clean_ion_global_map) {
  nrn_cleanup(clean_ion_global_map);
}

void DataLoader::LoadCoreneuronData(int argc, char **argv,
                                    bool nrnmpi_under_nrncontrol,
                                    bool run_setup_cleanup) {
  // nrnmpi_under_nrncontrol=true allows parallel data loading without "-m"
  // flag, see rnmpi_init()
  nrn_init_and_load_data(argc, argv, nrnmpi_under_nrncontrol,
                         run_setup_cleanup);
}

int DataLoader::GetMyNrnThreadsCount() {
  assert(nrn_threads);
  return std::accumulate(nrn_threads, nrn_threads + nrn_nthread, 0,
                         [](int n, NrnThread &nt) { return n + nt.ncell; });
}

hpx_action_t DataLoader::InitMechanisms = 0;
int DataLoader::InitMechanisms_handler() {
  NEUROX_MEM_PIN(uint64_t);

  // To insert mechanisms in the right order, we must first calculate
  // dependencies
  int my_nrn_threads_count = GetMyNrnThreadsCount();

  if (my_nrn_threads_count == 0) return neurox::wrappers::MemoryUnpin(target);

  for (int i = 0; i < my_nrn_threads_count; i++) {
    assert(nrn_threads[i].ncell == 1);
  }

  // Different nrn_threads[i] have diff mechanisms; we'll get the union of all
  // neurons' mechs
  std::list<NrnThreadMembList *> ordered_mechs;  // list of unique mechanisms
  std::set<int> unique_mech_ids;                 // list of unique mechanism ids

  // insert all mechs from first neuron
  for (NrnThreadMembList *tml = nrn_threads[0].tml; tml != NULL;
       tml = tml->next) {
    ordered_mechs.push_back(tml);
    unique_mech_ids.insert(tml->index);
  }

  // insert all mechs from other neurons that do not yet exist, in the right
  // order
  for (int i = 1; i < my_nrn_threads_count; i++)
    for (NrnThreadMembList *tml = nrn_threads[i].tml; tml->next != NULL;
         tml = tml->next)
      if (unique_mech_ids.find(tml->next->index) ==
          unique_mech_ids.end())  // if next mech does not exist
      {  // find correct position in list and insert it there:
        for (auto it = ordered_mechs.begin(); it != ordered_mechs.end();
             it++)  //...for all existing mechs
          if ((*it)->index == tml->index) {
            auto it_next = std::next(it, 1);
            ordered_mechs.insert(it_next, tml->next);  // reminder: .insert adds
                                                       // elements in position
                                                       // before iterator
            // we have it -> it_new -> it_next. Now we will set the value of
            // next pointer
            // auto it_new = std::prev(it_next,1);
            //(*it)->next = *it_new; //Do not change next pointers or neuron is
            // incorrect
            //(*it_new)->next = *it_next;
            unique_mech_ids.insert(tml->next->index);
            break;
          }
      }
  assert(unique_mech_ids.size() == ordered_mechs.size());

  std::vector<int> mech_ids_serial;
  std::vector<int> dependencies_count_serial;
  std::vector<int> successors_count_serial;
  std::vector<int> dependencies_serial;
  std::vector<int> successors_serial;

  for (auto tml_it = ordered_mechs.begin(); tml_it != ordered_mechs.end();
       tml_it++) {
    auto &tml = *tml_it;
    int type = tml->index;

    vector<int> successors;
    vector<int> dependencies;
    assert(nrn_watch_check[type] == NULL);  // not supported yet

    if (input_params_->mechs_parallelism_) {
      for (auto &tml2 : ordered_mechs) {
        int otherType = tml2->index;
        for (int d = 0; d < nrn_prop_dparam_size_[otherType]; d++) {
          int ptype = memb_func[otherType].dparam_semantics[d];
          if (otherType == type && ptype > 0 && ptype < 1000) {
            if (std::find(dependencies.begin(), dependencies.end(), ptype) ==
                dependencies.end())
              dependencies.push_back(ptype);  // parent on dependency graph
          }
          if (otherType != type && ptype == type) {
            if (std::find(successors.begin(), successors.end(), otherType) ==
                successors.end())
              successors.push_back(otherType);  // children on dependency graph
          }
        }
      }
    } else {
      // all except last have one successor
      if (tml->index != ordered_mechs.back()->index) {
        auto tml_next_it = std::next(tml_it, 1);
        successors.push_back((*tml_next_it)->index);
      }
      // all except first have one dependency
      if (tml->index != CAP) {
        auto tml_prev_it = std::prev(tml_it, 1);
        dependencies.push_back((*tml_prev_it)->index);
      }
    }

    // serialize data
    mech_ids_serial.push_back(type);
    dependencies_count_serial.push_back(dependencies.size());
    successors_count_serial.push_back(successors.size());
    dependencies_serial.insert(dependencies_serial.end(), dependencies.begin(),
                               dependencies.end());
    successors_serial.insert(successors_serial.end(), successors.begin(),
                             successors.end());
  }

  if (input_params_->pattern_stim_[0] != '\0')  //"initialized"
  {
    /*not an error: should be initialized already by coreneuron
     * (and above). in the future this should be contidtional
     * (once we get rid of coreneuron data loading) */
    assert(0);
  }

  // set mechanisms dependencies
  if (neurox::ParallelExecution() &&
      input_params_->branch_parallelism_complexity_ > 0) {
    /* broadcast dependencies, most complete dependency graph will be used
     * across the network (this solves issue of localities loading morphologies
     * without all mechanisms, and processing branches of other localities where
     * those missing mechanisms exist) */
    hpx_bcast_rsync(
        DataLoader::SetMechanisms, mech_ids_serial.data(),
        sizeof(int) * mech_ids_serial.size(), dependencies_count_serial.data(),
        sizeof(int) * dependencies_count_serial.size(),
        dependencies_serial.data(), sizeof(int) * dependencies_serial.size(),
        successors_count_serial.data(),
        sizeof(int) * successors_count_serial.size(), successors_serial.data(),
        sizeof(int) * successors_serial.size());
  } else {
    // regular setting of mechanisms: all localities have the dependency graph
    // for their morphologies
    DataLoader::SetMechanisms2(
        mech_ids_serial.size(), mech_ids_serial.data(),
        dependencies_count_serial.data(), dependencies_serial.data(),
        successors_count_serial.data(), successors_serial.data());
  }
  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t DataLoader::Init = 0;
int DataLoader::Init_handler() {
  NEUROX_MEM_PIN(uint64_t);
  all_neurons_gids_ = new std::vector<int>();
  locality_mutex_ = hpx_lco_sema_new(1);

  // even without load balancing, we may require the benchmark info for
  // outputting statistics
  if (hpx_get_my_rank() == 0 &&
      (input_params_->load_balancing_ || input_params_->output_statistics_))
    load_balancing_ = new tools::LoadBalancing();

  if (neurox::ParallelExecution()  // disable output of netcons for parallel
                                   // loading
      && input_params_->output_netcons_dot) {
    input_params_->output_netcons_dot = false;
    if (hpx_get_my_rank() == 0)
      printf("Warning: output of netcons.dot disabled for parallel loading\n");
  }

  if (input_params_->output_netcons_dot) {
    assert(HPX_LOCALITIES == 1);
    file_netcons_ = fopen(string("netcons.dot").c_str(), "wt");
    fprintf(file_netcons_, "digraph G\n{ bgcolor=%s; layout=circo;\n",
            "transparent");
#if NEUROX_INPUT_DATALOADER_OUTPUT_EXTERNAL_NETCONS == true
    fprintf(file_netcons_, "external [color=gray fontcolor=gray];\n");
#endif
  }

  // initiate map of locality to branch netcons (if needed)
  if (input_params_->locality_comm_reduce_) {
    assert(locality::netcons_branches_ == nullptr);
    locality::neurons_ = new vector<hpx_t>();
    locality::netcons_branches_ = new map<neuron_id_t, vector<hpx_t>>();
    locality::netcons_somas_ = new map<neuron_id_t, vector<hpx_t>>();
  }

  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t DataLoader::InitNeurons = 0;
int DataLoader::InitNeurons_handler() {
  NEUROX_MEM_PIN(uint64_t);

  int my_nrn_threads_count = GetMyNrnThreadsCount();

  if (my_nrn_threads_count == 0) return neurox::wrappers::MemoryUnpin(target);

#if NEUROX_INPUT_DATALOADER_OUTPUT_CORENEURON_COMPARTMENTS == true
  if (input_params_->output_compartments_dot_) {
    for (int i = 0; i < my_nrn_threads_count; i++) {
      neuron_id_t neuron_id = GetNeuronIdFromNrnThreadId(i);
      FILE *file_compartments = fopen(
          string("compartments_" + to_string(neuron_id) + "_NrnThread.dot")
              .c_str(),
          "wt");
      fprintf(file_compartments, "graph G%d\n{  node [shape=cylinder];\n",
              neuron_id);

      // for all nodes in this NrnThread
      NrnThread *nt = &nrn_threads[i];
      for (int n = nt->ncell; n < nt->end; n++)
        fprintf(file_compartments, "%d -- %d;\n", nt->_v_parent_index[n], n);
      fprintf(file_compartments, "}\n");
      fclose(file_compartments);
    }
  }
#endif

  my_neurons_addr_ = new std::vector<hpx_t>();
  my_neurons_gids_ = new std::vector<int>();

  hpx_par_for_sync(DataLoader::CreateNeuron, 0, my_nrn_threads_count, nullptr);

  // all neurons created, advertise them
  int my_rank = hpx_get_my_rank();
  hpx_bcast_rsync(
      DataLoader::AddNeurons, my_neurons_gids_->data(),
      sizeof(int) * my_neurons_gids_->size(), my_neurons_addr_->data(),
      sizeof(hpx_t) * my_neurons_addr_->size(), &my_rank, sizeof(my_rank));

  my_neurons_gids_->clear();
  delete my_neurons_gids_;
  my_neurons_gids_ = nullptr;
  my_neurons_addr_->clear();
  delete my_neurons_addr_;
  my_neurons_addr_ = nullptr;

  return neurox::wrappers::MemoryUnpin(target);
}

hpx_action_t DataLoader::AddNeurons = 0;
int DataLoader::AddNeurons_handler(const int nargs, const void *args[],
                                   const size_t sizes[]) {
  /**
   * nargs=3 where
   * args[0] = neurons Gids
   * args[1] = neurons hpx addr
   * args[2] = sender rank
   */

  NEUROX_MEM_PIN(uint64_t);
  assert(nargs == 3);
  assert(sizes[0] / sizeof(int) == sizes[1] / sizeof(hpx_t));

  const int recv_neurons_count = sizes[0] / sizeof(int);
  const int *neurons_gids = (const int *)args[0];
  const hpx_t *neurons_addr = (const hpx_t *)args[1];

  hpx_lco_sema_p(locality_mutex_);

  all_neurons_gids_->insert(all_neurons_gids_->end(), neurons_gids,
                            neurons_gids + recv_neurons_count);

  // add the values to neurox::neurons
  hpx_t *neurons_new = new hpx_t[neurox::neurons_count_ + recv_neurons_count];
  copy(neurox::neurons_, neurox::neurons_ + neurox::neurons_count_,
       neurons_new);
  copy(neurons_addr, neurons_addr + recv_neurons_count,
       neurons_new + neurox::neurons_count_);

  delete[] neurox::neurons_;
  neurox::neurons_count_ += recv_neurons_count;
  neurox::neurons_ = neurons_new;

  assert(all_neurons_gids_->size() == neurox::neurons_count_);

  // initialize local neurons
  const int sender_rank = *(const int *)args[2];

  if (sender_rank == hpx_get_my_rank())  // if these are my neurons
  {
    if (input_params_->locality_comm_reduce_) {
      assert(locality::neurons_->size() == 0);
      locality::neurons_->insert(locality::neurons_->end(), neurons_addr,
                                 neurons_addr + recv_neurons_count);
      locality::neurons_->shrink_to_fit();
    }
  }
  hpx_lco_sema_v_sync(locality_mutex_);
  NEUROX_MEM_UNPIN;
}

hpx_action_t DataLoader::SetMechanisms = 0;
int DataLoader::SetMechanisms_handler(const int nargs, const void *args[],
                                      const size_t sizes[]) {
  /**
   * nargs=5 where
   * args[0] = mechanisms ids (in the correct dependency order)
   * args[1] = successors count per mechanism
   * args[2] = successors ids
   * args[3] = dependendencies count per mechanism
   * args[4] = dependendencies ids
   */

  NEUROX_MEM_PIN(uint64_t);
  assert(nargs == 5);

  const int *ordered_mechs = (const int *)args[0];
  const int *dependencies_count = (const int *)args[1];
  const int *dependencies = (const int *)args[2];
  const int *successors_count = (const int *)args[3];
  const int *successors = (const int *)args[4];

  const int mechs_count = sizes[0] / sizeof(int);

  SetMechanisms2(mechs_count, ordered_mechs, dependencies_count, dependencies,
                 successors_count, successors);

  return neurox::wrappers::MemoryUnpin(target);
}

void DataLoader::SetMechanisms2(const int mechs_count, const int *mech_ids,
                                const int *dependencies_count,
                                const int *dependencies,
                                const int *successors_count,
                                const int *successors) {
  hpx_lco_sema_p(locality_mutex_);

  // if I'm receiving new information, rebuild mechanisms
  if (mechs_count > neurox::mechanisms_count_) {
    // delete existing data (if any)
    if (neurox::mechanisms_count_ > 0) {
      for (int m = 0; m < neurox::mechanisms_count_; m++) delete mechanisms_[m];
      delete[] neurox::mechanisms_;
      delete[] neurox::mechanisms_map_;
    }

    neurox::mechanisms_count_ = mechs_count;
    neurox::mechanisms_map_ = new int[n_memb_func];
    neurox::mechanisms_ = new Mechanism *[mechs_count];

    for (int i = 0; i < n_memb_func; i++) neurox::mechanisms_map_[i] = -1;

    int dependencies_offset = 0;
    int successors_offset = 0;
    for (int m = 0; m < mechs_count; m++) {
      int type = mech_ids[m];
      neurox::mechanisms_map_[type] = m;

      int sym_length =
          memb_func[type].sym ? std::strlen(memb_func[type].sym) : 0;
      const int *mech_dependencies = dependencies_count[m] == 0
                                         ? nullptr
                                         : &dependencies[dependencies_offset];
      const int *mech_successors =
          successors_count[m] == 0 ? nullptr : &successors[successors_offset];
      Mechanism *mech = new Mechanism(
          type, nrn_prop_param_size_[type], nrn_prop_dparam_size_[type],
          nrn_is_artificial_[type], pnt_map[type], nrn_is_ion(type), sym_length,
          memb_func[type].sym, memb_func[type], dependencies_count[m],
          mech_dependencies, successors_count[m], mech_successors);
      neurox::mechanisms_[m] = mech;

      successors_offset += successors_count[m];
      dependencies_offset += dependencies_count[m];
    }

    // set parent ion index
    for (int m = 0; m < mechs_count; m++) {
      if (input_params_->mechs_parallelism_) {
        Mechanism *mech = neurox::mechanisms_[m];
        for (int d = 0; d < mech->dependencies_count_; d++) {
          Mechanism *parent = GetMechanismFromType(mech->dependencies_[d]);
          if (strcmp("SK_E2", mech->memb_func_.sym) == 0 &&
              strcmp("ca_ion", parent->memb_func_.sym) == 0)
            continue;  // TODO hard coded exception

          if (parent->GetIonIndex() < Mechanism::IonTypes::kSizeWriteableIons)
            mech->dependency_ion_index_ = parent->GetIonIndex();
        }
      }
    }
  }
  hpx_lco_sema_v_sync(locality_mutex_);
}

hpx_action_t DataLoader::Finalize = 0;
int DataLoader::Finalize_handler() {
  NEUROX_MEM_PIN(uint64_t);

  if (input_params_->output_netcons_dot) {
    fprintf(file_netcons_, "}\n");
    fclose(file_netcons_);
  }

  if (input_params_->output_mechanisms_dot_) {
    FILE *file_mechs =
        fopen(string("mechanisms_" + std::to_string(hpx_get_my_rank()) + ".dot")
                  .c_str(),
              "wt");
    fprintf(
        file_mechs, "digraph G\n{ bgcolor=%s; %s\n", "transparent",
        !input_params_->mechs_parallelism_ ? "layout=circo; scale=0.23;" : "");
    fprintf(file_mechs, "graph [ratio=0.3];\n");
    fprintf(file_mechs, "%s [style=filled, shape=Mdiamond, fillcolor=beige];\n",
            "start");
    fprintf(file_mechs, "%s [style=filled, shape=Mdiamond, fillcolor=beige];\n",
            "end");
    fprintf(file_mechs, "\"%s (%d)\" [style=filled, fillcolor=beige];\n",
            GetMechanismFromType(CAP)->memb_func_.sym, CAP);
    if (!input_params_->mechs_parallelism_) {
      fprintf(file_mechs, "end -> start [color=transparent];\n");
      fprintf(file_mechs, "start -> \"%s (%d)\";\n",
              GetMechanismFromType(CAP)->memb_func_.sym, CAP);
    }
    for (int m = 0; m < mechanisms_count_; m++) {
      Mechanism *mech = neurox::mechanisms_[m];

      if (mech->pnt_map_ > 0)  // if is point process make it dotted
        fprintf(file_mechs, "\"%s (%d)\" [style=dashed];\n",
                mech->memb_func_.sym, mech->type_);

      if (mech->dependencies_count_ == 0 &&
          mech->type_ != CAP)  // top mechanism
        fprintf(file_mechs, "%s -> \"%s (%d)\";\n", "start",
                mech->memb_func_.sym, mech->type_);

      if (mech->successors_count_ == 0 &&
          mech->type_ != CAP)  // bottom mechanism
        fprintf(file_mechs, "\"%s (%d)\" -> %s;\n", mech->memb_func_.sym,
                mech->type_, "end");

      for (int s = 0; s < mech->successors_count_; s++) {
        Mechanism *successor = GetMechanismFromType(mech->successors_[s]);
        fprintf(file_mechs, "\"%s (%d)\" -> \"%s (%d)\";\n",
                mech->memb_func_.sym, mech->type_, successor->memb_func_.sym,
                successor->type_);
      }

      if (input_params_->mechs_parallelism_)
        for (int d = 0; d < mech->dependencies_count_; d++) {
          Mechanism *parent = GetMechanismFromType(mech->dependencies_[d]);
          if (strcmp("SK_E2", mech->memb_func_.sym) == 0 &&
              strcmp("ca_ion", parent->memb_func_.sym) == 0)
            continue;  // TODO: hardcoded exception
          if (parent->GetIonIndex() <
              Mechanism::IonTypes::kSizeWriteableIons)  // ie is writeable
            fprintf(
                file_mechs,
                "\"%s (%d)\" -> \"%s (%d)\" [style=dashed, arrowtype=open];\n",
                mech->memb_func_.sym, mech->type_, parent->memb_func_.sym,
                parent->type_);
        }
    }
    fprintf(file_mechs, "}\n");
    fclose(file_mechs);
  }

#ifndef NDEBUG
  if (HPX_LOCALITY_ID == 0) {
    for (int m = 0; m < neurox::mechanisms_count_; m++) {
      Mechanism *mech = neurox::mechanisms_[m];
      printf(
          "- %s (%d), dataSize %d, pdataSize %d, isArtificial %d, pntMap %d, "
          "isIon %d, symLength %d, %d successors, %d dependencies, %d state "
          "vars\n",
          mech->memb_func_.sym, mech->type_, mech->data_size_,
          mech->pdata_size_, mech->is_artificial_, mech->pnt_map_,
          mech->is_ion_, mech->sym_length_, mech->successors_count_,
          mech->dependencies_count_,
          mech->state_vars_ ? mech->state_vars_->count_ : 0);
    }
  }
#endif

  all_neurons_gids_->clear();
  delete all_neurons_gids_;
  all_neurons_gids_ = nullptr;

  nrn_setup_cleanup();
  hpx_lco_delete_sync(locality_mutex_);

#if defined(NDEBUG)
  // if not on debug, there's no CoreNeuron comparison, so data can be
  // cleaned-up now, except global_ion_map
  DataLoader::CleanCoreneuronData(false);
#endif

  if (input_params_->output_statistics_)
    tools::LoadBalancing::LoadBalancing::PrintTable();

  delete load_balancing_;
  load_balancing_ = nullptr;

  // delete duplicates in locality map of netcons to addr
  if (input_params_->locality_comm_reduce_ && locality::netcons_branches_) {
    for (auto &map_it : (*locality::netcons_branches_)) {
      vector<hpx_t> &addrs = map_it.second;
      std::sort(addrs.begin(), addrs.end());
      addrs.erase(unique(addrs.begin(), addrs.end()), addrs.end());
    }

    for (auto &map_it : (*locality::netcons_somas_)) {
      vector<hpx_t> &addrs = map_it.second;
      std::sort(addrs.begin(), addrs.end());
      addrs.erase(unique(addrs.begin(), addrs.end()), addrs.end());
    }
  }
  return neurox::wrappers::MemoryUnpin(target);
}

void DataLoader::GetSubSectionFromCompartment(deque<Compartment *> &sub_section,
                                              Compartment *top_compartment) {
  // parent added first, to respect solvers logic, of parents' ids first
  sub_section.push_back(top_compartment);
  for (int c = 0; c < top_compartment->branches_.size(); c++)
    GetSubSectionFromCompartment(sub_section, top_compartment->branches_.at(c));
}

void DataLoader::GetMechInstanceMap(
    const deque<Compartment *> &compartments,
    vector<map<int, int>> &mech_instances_map) {
  vector<deque<int>> mech_instances_ids(
      neurox::mechanisms_count_);  // mech offset -> list of mech instance id
  for (const Compartment *comp : compartments)
    for (int m = 0; m < comp->mechs_types_.size(); m++)  // for all instances
    {
      int type = comp->mechs_types_.at(m);
      int mech_offset = neurox::mechanisms_map_[type];
      mech_instances_ids[mech_offset].push_back(comp->mechs_instances_[m]);
    }

  // convert neuron mech-instances ids from neuron- to branch-level
  for (int m = 0; m < neurox::mechanisms_count_; m++)
    for (int i = 0; i < mech_instances_ids[m].size(); i++) {
      int oldInstanceId = mech_instances_ids.at(m).at(i);
      mech_instances_map[m][oldInstanceId] = i;
    }
}

void DataLoader::GetNetConsBranchData(
    const deque<Compartment *> &compartments, vector<NetconX> &branch_netcons,
    vector<neuron_id_t> &branch_netcons_pre_id,
    vector<floble_t> &branch_weights,
    vector<map<int, int>> *mech_instances_map) {
  branch_netcons.clear();
  branch_netcons_pre_id.clear();
  branch_weights.clear();

  for (const auto &comp : compartments) {
    branch_netcons.insert(branch_netcons.end(), comp->netcons_.begin(),
                          comp->netcons_.end());
    branch_netcons_pre_id.insert(branch_netcons_pre_id.end(),
                                 comp->netcons_pre_syn_ids_.begin(),
                                 comp->netcons_pre_syn_ids_.end());
    branch_weights.insert(branch_weights.end(), comp->netcons_weights_.begin(),
                          comp->netcons_weights_.end());
  }

  // correct weighIndex variables
  int weight_offset = 0;
  for (NetconX &nc : branch_netcons) {
    nc.weight_index_ = weight_offset;
    weight_offset += nc.weights_count_;
  }

  // convert mech instance id from neuron to branch level
  if (mech_instances_map)
    for (NetconX &nc : branch_netcons)
      nc.mech_instance_ =
          (*mech_instances_map)[neurox::mechanisms_map_[nc.mech_type_]]
                              [nc.mech_instance_];
}

void DataLoader::GetVecPlayBranchData(
    const deque<Compartment *> &compartments, vector<floble_t> &vecplay_t_data,
    vector<floble_t> &vecplay_y_data, vector<PointProcInfo> &vecplay_info,
    vector<map<int, int>> *mech_instances_map) {
  vecplay_t_data.clear();
  vecplay_y_data.clear();
  vecplay_info.clear();

  // convert node id and mech instance id in PointProcess from neuron to branch
  // level
  if (mech_instances_map) {
    std::map<int, int> from_old_to_new_compartment_id;
    for (int n = 0; n < compartments.size(); n++)
      from_old_to_new_compartment_id[compartments.at(n)->id_] = n;

    for (int p = 0; p < vecplay_info.size(); p++) {
      PointProcInfo &ppi = vecplay_info[p];
      ppi.mech_instance =
          (offset_t)(*mech_instances_map)[neurox::mechanisms_map_[ppi.mech_type]]
                                        [ppi.mech_instance];
      ppi.node_id = from_old_to_new_compartment_id[ppi.node_id];
    }
  }

  for (const auto comp : compartments) {
    vecplay_t_data.insert(vecplay_t_data.end(), comp->vecplay_tdata_.begin(),
                          comp->vecplay_tdata_.end());
    vecplay_y_data.insert(vecplay_y_data.end(), comp->vecplay_ydata_.begin(),
                          comp->vecplay_ydata_.end());
    vecplay_info.insert(vecplay_info.end(), comp->vecplay_info_.begin(),
                        comp->vecplay_info_.end());
  }
}

int DataLoader::GetBranchData(
    const deque<Compartment *> &compartments, vector<floble_t> &data,
    vector<offset_t> &pdata, vector<unsigned char> &vdata, vector<offset_t> &p,
    vector<offset_t> &instances_count, vector<offset_t> &nodes_indices, int N,
    vector<DataLoader::IonInstancesInfo> &ions_instances_info,
    vector<map<int, int>> *mech_instance_map) {
  for (const auto comp : compartments) {
    assert(comp != NULL);
  }

  data.clear();
  pdata.clear();
  vdata.clear();
  p.clear();
  nodes_indices.clear();

  for (int m = 0; m < neurox::mechanisms_count_; m++) instances_count.at(m) = 0;

  int n = 0;  // number of compartments
  int vdata_pointer_offset = 0;

  ////// Basic information for RHS, D, A, B, V and area
  for (const auto comp : compartments) data.push_back(comp->rhs_);
  for (const auto comp : compartments) data.push_back(comp->d_);
  for (const auto comp : compartments) data.push_back(comp->a_);
  for (const auto comp : compartments) data.push_back(comp->b_);
  for (const auto comp : compartments) data.push_back(comp->v_);
  for (const auto comp : compartments) data.push_back(comp->area_);
  for (const auto comp : compartments) p.push_back(comp->p_);

  ////// Tree of neurons: convert from neuron- to branch-level
  std::map<int, int> from_old_to_new_compartment_id;
  for (const Compartment *comp : compartments)
    from_old_to_new_compartment_id[comp->id_] = n++;

  if (mech_instance_map) {
    p.at(0) = 0;  // top node gets parent Id 0 as in Coreneuron
    for (int i = 1; i < p.size(); i++)
      p.at(i) = from_old_to_new_compartment_id.at(p.at(i));
  }

  ////// Mechanisms instances: merge all instances of all compartments into
  /// instances of the branch
  vector<vector<floble_t>> data_mechs(neurox::mechanisms_count_);
  vector<vector<offset_t>> pdata_mechs(neurox::mechanisms_count_);
  vector<vector<offset_t>> nodes_indices_mechs(neurox::mechanisms_count_);
  vector<vector<unsigned char>> vdata_mechs(neurox::mechanisms_count_);
  vector<deque<int>> mechs_instances_ids(
      neurox::mechanisms_count_);  // mech-offset -> list of pointers to mech
                                   // instance value

  // from pair of <ion mech type, OLD node id> to ion offset in NEW
  // representation
  map<pair<int, offset_t>, offset_t> ion_instance_to_data_offset;

  for (const Compartment *comp : compartments) {
    int comp_data_offset = 0;
    int comp_pdata_offset = 0;
    int comp_vdata_offset = 0;
    for (int m = 0; m < comp->mechs_types_.size(); m++)  // for all instances
    {
      int type = comp->mechs_types_[m];
      int mech_offset = mechanisms_map_[type];
      assert(mech_offset >= 0 && mech_offset < mechanisms_count_);
      Mechanism *mech = mechanisms_[mech_offset];
      data_mechs[mech_offset].insert(
          data_mechs[mech_offset].end(), &comp->data[comp_data_offset],
          &comp->data[comp_data_offset + mech->data_size_]);
      pdata_mechs[mech_offset].insert(
          pdata_mechs[mech_offset].end(), &comp->pdata[comp_pdata_offset],
          &comp->pdata[comp_pdata_offset + mech->pdata_size_]);
      int new_id = from_old_to_new_compartment_id.at(comp->id_);
      nodes_indices_mechs[mech_offset].push_back(new_id);
      instances_count[mech_offset]++;
      comp_data_offset += mech->data_size_;
      comp_pdata_offset += mech->pdata_size_;
      if (mech->pnt_map_ > 0 || mech->vdata_size_ > 0) {
        assert((type == MechanismTypes::kIClamp && mech->vdata_size_ == 1 &&
                mech->pdata_size_ == 2 && mech->pnt_map_ > 0) ||
               (type == MechanismTypes::kStochKv && mech->vdata_size_ == 1 &&
                mech->pdata_size_ == 5 && mech->pnt_map_ == 0) ||
               ((type == MechanismTypes::kProbAMPANMDA_EMS ||
                 type == MechanismTypes::kProbGABAAB_EMS) &&
                mech->vdata_size_ == 2 && mech->pdata_size_ == 3 &&
                mech->pnt_map_ > 0));

        size_t total_vdata_size = 0;
        if (type == MechanismTypes::kIClamp ||
            type == MechanismTypes::kProbAMPANMDA_EMS ||
            type == MechanismTypes::kProbGABAAB_EMS)
          total_vdata_size += sizeof(Point_process);
        if (type == MechanismTypes::kStochKv ||
            type == MechanismTypes::kProbAMPANMDA_EMS ||
            type == MechanismTypes::kProbGABAAB_EMS)
          total_vdata_size += sizeof(nrnran123_State);
        vdata_mechs[mech_offset].insert(
            vdata_mechs[mech_offset].end(), &comp->vdata_[comp_vdata_offset],
            &comp->vdata_[comp_vdata_offset + total_vdata_size]);
        comp_vdata_offset += total_vdata_size;
      }
    }
  }

  // merge all mechanisms vectors in the final one
  // store the offset of each mechanism data (for later)
  for (int m = 0; m < neurox::mechanisms_count_; m++) {
    Mechanism *mech = neurox::mechanisms_[m];
    int data_offset = 0;
    int pdata_offset = 0;
    int vdata_offset = 0;
    for (int i = 0; i < nodes_indices_mechs[m].size();
         i++)  // for all instances
    {
      assert(ion_instance_to_data_offset.find(
                 make_pair(mech->type_, nodes_indices_mechs[m][i])) ==
             ion_instance_to_data_offset.end());
      if (mech->is_ion_ &&
          input_params_->branch_parallelism_complexity_ >
              0)  // for pdata calculation
        ion_instance_to_data_offset[make_pair(
            mech->type_, nodes_indices_mechs[m][i])] = data.size();
      data.insert(data.end(), &data_mechs[m][data_offset],
                  &data_mechs[m][data_offset + mech->data_size_]);
      nodes_indices.push_back(nodes_indices_mechs[m][i]);
      data_offset += mech->data_size_;

      if (mech->pnt_map_ > 0 || mech->vdata_size_ > 0) {
        int total_vdata_size = 0;
        if (mech->type_ == MechanismTypes::kIClamp ||
            mech->type_ == MechanismTypes::kProbAMPANMDA_EMS ||
            mech->type_ == MechanismTypes::kProbGABAAB_EMS)
          total_vdata_size += sizeof(Point_process);
        if (mech->type_ == MechanismTypes::kStochKv ||
            mech->type_ == MechanismTypes::kProbAMPANMDA_EMS ||
            mech->type_ == MechanismTypes::kProbGABAAB_EMS)
          total_vdata_size += sizeof(nrnran123_State);
        assert(total_vdata_size > 0);
        vdata.insert(vdata.end(), &vdata_mechs[m][vdata_offset],
                     &vdata_mechs[m][vdata_offset + total_vdata_size]);
        vdata_offset += total_vdata_size;
      }

      if (input_params_->branch_parallelism_complexity_ >
          0)  // if we need to recalculate offsets or remove padding
      {
        for (int p = pdata_offset; p < pdata_offset + mech->pdata_size_; p++) {
          offset_t pd = pdata_mechs.at(m).at(p);
          int ptype = memb_func[mech->type_].dparam_semantics[p - pdata_offset];
          switch (ptype) {
            case -1:  //"area" (6th field)
            {
              assert(pd >= N * 5 && pd < N * 6);
              offset_t old_id = pd - N * 5;
              offset_t new_id = from_old_to_new_compartment_id.at(old_id);
              pdata_mechs.at(m).at(p) = n * 5 + new_id;
              break;
            }
            case -2:  //"iontype"
              // do nothing, its a flag (the 'iontype', see nrnoc/eion.c)
              break;
            case -3:      //"cvodeieq"
            case -5:      //"pointer"
              assert(0);  // not used
              break;
            case -4:  //"netsend"
            case -6:  //"pntproc"
            case -7:  //"bbcorepointer"
              pdata_mechs.at(m).at(p) = (offset_t)vdata_pointer_offset++;
              break;
            case -8:  // watch condition, not supported
              assert(0);
              break;
            default:
              if (ptype > 0 && ptype < 1000)  // name preffixed by '#'
              {
                // ptype is the ion (mechanism) type it depends on
                // pdata is an offset in nt->data (a var in the ion)
                // pd points to SoA notation, independently of the LAYOUT
                // (converted before)

                Mechanism *ion = neurox::GetMechanismFromType(ptype);
                Mechanism::IonTypes ion_offset = ion->GetIonIndex();
                IonInstancesInfo &ion_info =
                    ions_instances_info.at((int)ion_offset);
                int data_start = ion_info.data_start;
                assert(pd >= data_start && pd < ion_info.data_end);

                int instance_offset =
                    trunc((double)(pd - data_start) / (double)ion->data_size_);
                int instance_variable_offset =
                    (pd - data_start) % ion->data_size_;
                int node_id = ion_info.node_ids.at(instance_offset);
                int new_node_id = from_old_to_new_compartment_id.at(node_id);
                pdata_mechs.at(m).at(p) =
                    ion_instance_to_data_offset.at(
                        make_pair(ion->type_, new_node_id)) +
                    instance_variable_offset;
                assert(pdata_mechs.at(m).at(p) >= n * 6);
              } else if (ptype >= 1000)  // name not preffixed
              {
                // (concentration = ptype-1000;) //do nothing: value of
                // concentration summed with 1000
              } else
                throw std::runtime_error("Unknown pdata type %d (FLAG3)\n");
              break;
          }
        }
      }
      pdata.insert(pdata.end(), &pdata_mechs[m][pdata_offset],
                   &pdata_mechs[m][pdata_offset + mech->pdata_size_]);
      pdata_offset += mech->pdata_size_;
    }
    data_mechs[m].clear();
    pdata_mechs[m].clear();
    vdata_mechs[m].clear();
    nodes_indices_mechs[m].clear();

    // convert neuron mech-instances ids from neuron- to branch-level
    if (mech_instance_map)
      for (int i = 0; i < mechs_instances_ids[m].size(); i++) {
        int old_instance_id = mechs_instances_ids[m][i];
        (*mech_instance_map)[m][old_instance_id] = i;
      }
  }
  return n;
}

bool CompareCompartmentPtrId(Compartment *a, Compartment *b) {
  return a->id_ < b->id_;
}

hpx_t DataLoader::CreateBranch(
    const int nrn_threadId, hpx_t soma_branch_addr,
    const deque<Compartment *> &all_compartments, Compartment *top_compartment,
    vector<DataLoader::IonInstancesInfo> &ions_instances_info,
    double max_work_per_section, int thvar_index /*AIS*/,
    floble_t ap_threshold /*AIS*/) {

  assert(top_compartment != NULL);
  offset_t n;              // number of compartments in branch
  vector<floble_t> data;   // compartments info (RHS, D, A, B, V, AREA)*n
  vector<offset_t> pdata;  // pointers to data
  vector<offset_t> p;      // parent nodes index
  vector<offset_t> instances_count(mechanisms_count_);
  vector<offset_t> nodes_indices;
  vector<hpx_t> branches;
  vector<unsigned char> vdata;  // Serialized Point Processes and Random123
  int N = all_compartments.size();

  bool is_soma = soma_branch_addr == HPX_NULL;

  // Vector Play instances
  vector<floble_t> vecplay_t;
  vector<floble_t> vecplay_y;
  vector<PointProcInfo> vecplay_info;

  // branch NetCons
  vector<NetconX> branch_netcons;
  vector<neuron_id_t> branch_netcons_pre_id;
  vector<floble_t> branch_weights;

  Compartment *bottom_compartment = nullptr;

  /* Benchmark this subsection, if necessary */
  double time_elapsed = -1;
  if (input_params_->load_balancing_ || input_params_->output_statistics_ ||
      input_params_->branch_parallelism_complexity_ > 0) {

    // allocate GAS memory for this temporary branch
    hpx_t temp_branch_addr =
      hpx_gas_alloc_local(1, sizeof(Branch), Vectorizer::kMemoryAlignment);
    bool run_benchmark_and_clear = true;
    int dumb_threshold_offset = 0;

    // Get execution time for the whole subsection
    deque<Compartment *> whole_sub_section;

    // mech-offset -> ( map[old instance]->to new instance )
    vector<map<int, int>> mech_instances_map(neurox::mechanisms_count_);
    GetMechInstanceMap(whole_sub_section, mech_instances_map);

    GetSubSectionFromCompartment(whole_sub_section, top_compartment);
    n = GetBranchData(whole_sub_section, data, pdata, vdata, p, instances_count,
                      nodes_indices, N, ions_instances_info, &mech_instances_map);
    GetVecPlayBranchData(whole_sub_section, vecplay_t, vecplay_y, vecplay_info,
                         &mech_instances_map);
    GetNetConsBranchData(whole_sub_section, branch_netcons, branch_netcons_pre_id,
                         branch_weights, &mech_instances_map);
    whole_sub_section.clear();

    // run benchmark //TODO disable mechanism parallelism
    hpx_call_sync(
        temp_branch_addr, Branch::Init, &time_elapsed,
        sizeof(time_elapsed),  // output
        &n, sizeof(offset_t), &nrn_threadId, sizeof(int),
        &dumb_threshold_offset, sizeof(int),
        data.size() > 0 ? data.data() : nullptr, sizeof(floble_t) * data.size(),
        pdata.size() > 0 ? pdata.data() : nullptr,
        sizeof(offset_t) * pdata.size(), instances_count.data(),
        instances_count.size() * sizeof(offset_t), nodes_indices.data(),
        nodes_indices.size() * sizeof(offset_t), &soma_branch_addr,
        sizeof(hpx_t), nullptr, 0,              // no branches
        p.data(), sizeof(offset_t) * p.size(),  // force use of parent index
        vecplay_t.size() > 0 ? vecplay_t.data() : nullptr,
        sizeof(floble_t) * vecplay_t.size(),
        vecplay_y.size() > 0 ? vecplay_y.data() : nullptr,
        sizeof(floble_t) * vecplay_y.size(),
        vecplay_info.size() > 0 ? vecplay_info.data() : nullptr,
        sizeof(PointProcInfo) * vecplay_info.size(),
        branch_netcons.size() > 0 ? branch_netcons.data() : nullptr,
        sizeof(NetconX) * branch_netcons.size(),
        branch_netcons_pre_id.size() > 0 ? branch_netcons_pre_id.data()
                                         : nullptr,
        sizeof(neuron_id_t) * branch_netcons_pre_id.size(),
        branch_weights.size() > 0 ? branch_weights.data() : nullptr,
        sizeof(floble_t) * branch_weights.size(),
        vdata.size() > 0 ? vdata.data() : nullptr,
        sizeof(unsigned char) * vdata.size(),
        &run_benchmark_and_clear, sizeof(bool));
    assert(time_elapsed > 0);
    hpx_gas_clear_affinity(temp_branch_addr);
  }

  if (input_params_->branch_parallelism_complexity_ == 0) {
    // Get information about the whole subsection
    n = GetBranchData(all_compartments, data, pdata, vdata, p, instances_count,
                      nodes_indices, N, ions_instances_info, NULL);
    GetVecPlayBranchData(all_compartments, vecplay_t, vecplay_y, vecplay_info,
                         NULL);
    GetNetConsBranchData(all_compartments, branch_netcons,
                         branch_netcons_pre_id, branch_weights, NULL);
  } else {
    deque<Compartment *> sub_section;

    // time_elapsed is the total execution time for this neuron
    // set max ammout of work (computation time) assigned to each subsection
    if (is_soma)
      max_work_per_section =
          time_elapsed / (wrappers::NumThreads() *
                          input_params_->branch_parallelism_complexity_);

    max_work_per_section-=1; //TODO HACK

    /*if this subsection does not exceed maximum time allowed per subsection*/
    if (time_elapsed <= max_work_per_section) {
      // subsection is the set of all children compartments (recursively)
      GetSubSectionFromCompartment(sub_section, top_compartment);
    } else {
      // otherwise, it's the top branch, plus all subregions
      sub_section.push_back(top_compartment);
      for (bottom_compartment = top_compartment;
           bottom_compartment->branches_.size() == 1;
           bottom_compartment = bottom_compartment->branches_.front())
        sub_section.push_back(bottom_compartment->branches_.front());
    }

    // mech-offset -> ( map[old instance]->to new instance )
    vector<map<int, int>> mech_instances_map(neurox::mechanisms_count_);
    GetMechInstanceMap(sub_section, mech_instances_map);

    n = GetBranchData(sub_section, data, pdata, vdata, p, instances_count,
                      nodes_indices, N, ions_instances_info,
                      &mech_instances_map);
    GetVecPlayBranchData(sub_section, vecplay_t, vecplay_y, vecplay_info,
                         &mech_instances_map);
    GetNetConsBranchData(sub_section, branch_netcons, branch_netcons_pre_id,
                         branch_weights, &mech_instances_map);
  }

  int neuron_rank = hpx_get_my_rank();
  if (input_params_->load_balancing_ || input_params_->output_statistics_) {
    if (input_params_->load_balancing_) {
      // ask master rank to query load balancing table and tell me where to
      // allocate this branch
      hpx_call_sync(HPX_THERE(0), tools::LoadBalancing::QueryLoadBalancingTable,
                    &neuron_rank, sizeof(int),       // output
                    &time_elapsed, sizeof(double));  // input [0]
    } else if (input_params_->output_statistics_)    // no load balancing
    {
      // compute node already decided, tell master rank to store benchmark info,
      // to be able to output statistics
      hpx_call_sync(HPX_THERE(0), tools::LoadBalancing::QueryLoadBalancingTable,
                    nullptr, 0,                     // output
                    &time_elapsed, sizeof(double),  // input[0]
                    &neuron_rank, sizeof(int));     // input[1]
    }

#ifndef NDEBUG
    printf("- %s %d, length %d, nrn_id %d, runtime %.6f ms, allocated to rank %d\n",
           is_soma ? "soma" : (thvar_index != -1 ? "AIS" : "dendrite"),
           top_compartment->id_, n, nrn_threadId, time_elapsed, neuron_rank);
#endif
  }

  // allocate and initialize branch on the respective owner
  hpx_t branch_addr = hpx_gas_alloc_local_at_sync(
      1, sizeof(Branch), Vectorizer::kMemoryAlignment, HPX_THERE(neuron_rank));

  // update hpx address of soma
  soma_branch_addr = is_soma ? branch_addr : soma_branch_addr;

  // allocate children branches recursively (if any)
  if (bottom_compartment) {
    for (size_t c = 0; c < bottom_compartment->branches_.size(); c++)
      branches.push_back(CreateBranch(
          nrn_threadId, soma_branch_addr, all_compartments,
          bottom_compartment->branches_[c], ions_instances_info,
          max_work_per_section, is_soma && c == 0 ? thvar_index - n : -1));
    /*offset in AIS = offset in soma - nt->end */

    if (is_soma) thvar_index = -1;
  }

  hpx_call_sync(
      branch_addr, Branch::Init, NULL, 0,  // no timing
      &n, sizeof(offset_t), &nrn_threadId, sizeof(int), &thvar_index,
      sizeof(int), data.size() > 0 ? data.data() : nullptr,
      sizeof(floble_t) * data.size(), pdata.size() > 0 ? pdata.data() : nullptr,
      sizeof(offset_t) * pdata.size(), instances_count.data(),
      instances_count.size() * sizeof(offset_t), nodes_indices.data(),
      nodes_indices.size() * sizeof(offset_t), &soma_branch_addr, sizeof(hpx_t),
      branches.size() ? branches.data() : nullptr,
      branches.size() ? sizeof(hpx_t) * branches.size() : 0,
      branches.size() ? nullptr : p.data(),
      branches.size() ? 0 : sizeof(offset_t) * p.size(),
      vecplay_t.size() > 0 ? vecplay_t.data() : nullptr,
      sizeof(floble_t) * vecplay_t.size(),
      vecplay_y.size() > 0 ? vecplay_y.data() : nullptr,
      sizeof(floble_t) * vecplay_y.size(),
      vecplay_info.size() > 0 ? vecplay_info.data() : nullptr,
      sizeof(PointProcInfo) * vecplay_info.size(),
      branch_netcons.size() > 0 ? branch_netcons.data() : nullptr,
      sizeof(NetconX) * branch_netcons.size(),
      branch_netcons_pre_id.size() > 0 ? branch_netcons_pre_id.data() : nullptr,
      sizeof(neuron_id_t) * branch_netcons_pre_id.size(),
      branch_weights.size() > 0 ? branch_weights.data() : nullptr,
      sizeof(floble_t) * branch_weights.size(),
      vdata.size() > 0 ? vdata.data() : nullptr,
      sizeof(unsigned char) * vdata.size());

  if (is_soma) { // create soma data
    int neuron_id = GetNeuronIdFromNrnThreadId(nrn_threadId);
    hpx_call_sync(branch_addr, Branch::InitSoma, NULL, 0, &neuron_id,
                  sizeof(neuron_id_t), &ap_threshold, sizeof(floble_t));

    hpx_lco_sema_p(locality_mutex_);
    my_neurons_addr_->push_back(branch_addr);
    my_neurons_gids_->push_back(neuron_id);
    hpx_lco_sema_v_sync(locality_mutex_);
  }
  return branch_addr;
}

hpx_action_t DataLoader::InitNetcons = 0;
int DataLoader::InitNetcons_handler() {
  NEUROX_MEM_PIN(Branch);
  NEUROX_RECURSIVE_BRANCH_ASYNC_CALL(DataLoader::InitNetcons);

  // TODO the fastest synaptic delay variable should be set here!
  if (local->soma_ && input_params_->output_netcons_dot)
    fprintf(file_netcons_, "%d [style=filled, shape=ellipse];\n",
            local->soma_->gid_);

  // set of <srcAddr, minDelay> synapses to notify
  std::deque<std::pair<hpx_t, floble_t>> netcons;

  // set of <src_gid, syn_min_delay> for dependencies
  std::deque<std::pair<int, spike_time_t>> dependencies;

  const floble_t impossibly_large_delay = 99999999;
  for (int i = 0; i < all_neurons_gids_->size(); i++) {
    neuron_id_t src_gid = all_neurons_gids_->at(i);

    // if I'm connected to it (ie is not artificial or non-existent)
    if (local->netcons_.find(src_gid) != local->netcons_.end()) {
      hpx_t src_addr = neurox::neurons_[i];
      floble_t min_delay = impossibly_large_delay;
      for (NetconX *&nc : local->netcons_.at(src_gid))
        if (nc->active_) min_delay = min(min_delay, nc->delay_);

#if NETCONS_OUTPUT_ADDITIONAL_VALIDATION_FILE == true
      if (input_params_->output_netcons_dot &&
          min_delay != impossibly_large_delay)
        fprintf(file_netcons_, "%d -> %d [label=\"%d (%.2fms)\"];\n", src_gid,
                local->soma->gid, local->netcons.at(src_gid).size(), min_delay);
#endif

      // tell the neuron to add a synapse to this branch and inform him of the
      // fastest netcon we have
      if (min_delay != impossibly_large_delay)  // if any active synapse
      {
        // add this netcon to list of netcons to be communicated
        netcons.push_back(make_pair(src_addr, min_delay));

        // add this pre-syn neuron as my time-dependency
        if (input_params_->synchronizer_ == SynchronizerIds::kBenchmarkAll ||
            input_params_->synchronizer_ == SynchronizerIds::kTimeDependency) {
          dependencies.push_back(make_pair(src_gid, min_delay));
        }
      }
    }
  }

  // inform pre-syn neuron that he connects to this branch or locality
  // (repeated entries will be removed by DataLoader::FilterLocalitySynapses)
  hpx_t netcons_lco = hpx_lco_and_new(netcons.size());
  int my_gid = local->soma_ ? local->soma_->gid_ : -1;
  hpx_t top_branch_addr =
      local->soma_ ? target : local->branch_tree_->top_branch_addr_;
  hpx_t top_branch_locality_addr =
      input_params_->locality_comm_reduce_ ? HPX_HERE : top_branch_addr;

  /* inform my pre-syn to add a synapse to the branch with address 'branch_addr'
   or to this locality addr (for locality-based reduction of communication) */
  hpx_t branch_addr = input_params_->locality_comm_reduce_ ? HPX_HERE : target;
  for (std::pair<hpx_t, floble_t> &nc : netcons)
    hpx_call(nc.first, DataLoader::AddSynapse, netcons_lco, &branch_addr,
             sizeof(hpx_t), &nc.second, sizeof(nc.second),
             &top_branch_locality_addr, sizeof(hpx_t), &my_gid, sizeof(int));
  hpx_lco_wait_reset(netcons_lco);
  hpx_lco_delete_sync(netcons_lco);

  /* inform my soma of my time dependencies. For locality-based reduction of
   * communication: when neuron spikes or steps, it will inform locality that
   * will then inform this soma */
  hpx_t dependencies_lco = hpx_lco_and_new(dependencies.size());
  bool init_phase = true;
  for (std::pair<int, spike_time_t> &dep : dependencies)
    hpx_call(top_branch_addr, TimeDependencySynchronizer::UpdateTimeDependency,
             dependencies_lco, &dep.first, sizeof(neuron_id_t), &dep.second,
             sizeof(spike_time_t), &init_phase, sizeof(bool));
  hpx_lco_wait_reset(dependencies_lco);
  hpx_lco_delete_sync(dependencies_lco);

  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT;
  NEUROX_MEM_UNPIN;
}

hpx_action_t DataLoader::AddSynapse = 0;
int DataLoader::AddSynapse_handler(const int nargs, const void *args[],
                                   const size_t[]) {
  NEUROX_MEM_PIN(Branch);
  assert(local->soma_);
  assert(nargs == 4);
  hpx_t addr = *(const hpx_t *)args[0];
  floble_t min_delay = *(const floble_t *)args[1];
  hpx_t soma_or_locality_addr = *(const hpx_t *)args[2];
  int destination_gid = *(const int *)args[3];
  Neuron::Synapse *syn = new Neuron::Synapse(
      addr, min_delay, soma_or_locality_addr, destination_gid);
  local->soma_->AddSynapse(syn);
  NEUROX_MEM_UNPIN
}

hpx_action_t DataLoader::FilterRepeatedLocalitySynapses = 0;
int DataLoader::FilterRepeatedLocalitySynapses_handler() {
  if (!input_params_->locality_comm_reduce_) return HPX_SUCCESS;

  NEUROX_MEM_PIN(Branch);
  /* this methods takes all outgoing synapses as <locality hpx_t, min_delay>,
   * removes all repeated localities (we have several copies), and uses only
   * the values of the min_delay per locality */

  // initial synapses, several per locality
  std::vector<Neuron::Synapse *> &synapses = local->soma_->synapses_;

  // filtered synapses (1 per locality)
  map<hpx_t, Neuron::Synapse *> synapses_loc;

  for (Neuron::Synapse *s : synapses) {
    if (synapses_loc.find(s->branch_addr_) == synapses_loc.end()) {
      synapses_loc[s->branch_addr_] = s;
      continue;
    }

    // synapse already exists, use only the one with min syn delay
    Neuron::Synapse *syn_old = synapses_loc.at(s->branch_addr_);
    assert(syn_old->branch_addr_ == s->branch_addr_);
    if (s->min_delay_ < syn_old->min_delay_)
      synapses_loc.at(s->branch_addr_) = s;
  }

  synapses.clear();
  for (auto &syn_it : synapses_loc) synapses.push_back(syn_it.second);
  assert(synapses.size() <= hpx_get_num_ranks());
  NEUROX_MEM_UNPIN
}

void DataLoader::RegisterHpxActions() {
  wrappers::RegisterZeroVarAction(
      DataLoader::FilterRepeatedLocalitySynapses,
      DataLoader::FilterRepeatedLocalitySynapses_handler);
  wrappers::RegisterZeroVarAction(DataLoader::Init, DataLoader::Init_handler);
  wrappers::RegisterZeroVarAction(DataLoader::InitMechanisms,
                                  DataLoader::InitMechanisms_handler);
  wrappers::RegisterZeroVarAction(DataLoader::InitNeurons,
                                  DataLoader::InitNeurons_handler);
  wrappers::RegisterZeroVarAction(DataLoader::InitNetcons,
                                  DataLoader::InitNetcons_handler);
  wrappers::RegisterZeroVarAction(DataLoader::Finalize,
                                  DataLoader::Finalize_handler);

  wrappers::RegisterMultipleVarAction(DataLoader::AddSynapse,
                                      DataLoader::AddSynapse_handler);
  wrappers::RegisterMultipleVarAction(DataLoader::AddNeurons,
                                      DataLoader::AddNeurons_handler);
  wrappers::RegisterMultipleVarAction(DataLoader::SetMechanisms,
                                      DataLoader::SetMechanisms_handler);
}
