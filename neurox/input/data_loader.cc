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
bool silence_p_ptr_warning = false;
bool silence_queue_item_warning = false;

// TODO: hard-coded procedures
int DataLoader::HardCodedVdataSize(int type) {
  assert(type != MechanismTypes::kInhPoissonStim);  // not supported yet
  int total_vdata_size = 0;

  if (HardCodedPntProcOffsetInPdata(type) != -1)
    total_vdata_size += sizeof(Point_process);
  if (HardCodedRNGOffsetInPdata(type) != -1)
    total_vdata_size += sizeof(nrnran123_State);
  if (HardCodedQueueItemOffsetInPdata(type) != -1)
    total_vdata_size += sizeof(void *);
  if (HardCodedPPtrOffsetInPdata(type) != -1)
    total_vdata_size += sizeof(void *);
  return total_vdata_size;
}

int DataLoader::HardCodedVdataCount(int type, char pnt_map) {
  assert(type != MechanismTypes::kInhPoissonStim);  // not supported yet
  int total_vdata_count = 0;
  if (pnt_map > 0) total_vdata_count++;
  if (HardCodedRNGOffsetInPdata(type) != -1) total_vdata_count++;
  if (HardCodedQueueItemOffsetInPdata(type) != -1) total_vdata_count++;
  if (HardCodedPPtrOffsetInPdata(type) != -1) total_vdata_count++;
  return total_vdata_count;
}

int DataLoader::HardCodedPntProcOffsetInPdata(int type) {
  // dictated by nrn_setup.cpp:
  // the area is always at ppvar[0]
  // the Point Process if always at position ppvar[1] --> pointing to vdata[0]
  if (GetMechanismFromType(type)->pnt_map_ > 0) return 1;
  return -1;  // means unexistent
}

int DataLoader::HardCodedPntProcOffsetInVdata(int type) {
  if (GetMechanismFromType(type)->pnt_map_ > 0)
    return 0;  // see comment in HardCodedPntProcOffsetInPdata
  return -1;   // means unexistent
}

int DataLoader::HardCodedRNGOffsetInPdata(int type) {
  assert(type != MechanismTypes::kInhPoissonStim);  // not supported yet

  if (type == MechanismTypes::kGluSynapse ||        //_p_rng_rel, ppvar[2]
      type == MechanismTypes::kNetStim ||           //_p_donotuse, ppvar[2]
      type == MechanismTypes::kProbAMPANMDA_EMS ||  //_p_rng_rel, ppvar[2]
      type == MechanismTypes::kProbGABAAB_EMS)      //_p_rng_rel, ppvar[2]
    return 2;
  if (type == MechanismTypes::kStochKv ||  //_p_rng, ppvar[3]
      type == MechanismTypes::kStochKv3)   //_p_rng, ppvar[3]
    return 3;
  return -1;  // means unexistent
}

int DataLoader::HardCodedRNGOffsetInVdata(int type) {
  assert(type != MechanismTypes::kInhPoissonStim);  // not supported yet

  int offset_rng_in_vdata = 0;
  // return type == MechanismTypes::kStochKv ? 0 : 1;
  if (GetMechanismFromType(type)->pnt_map_ > 0)
    offset_rng_in_vdata++;  // include offset for Point_process*
  return offset_rng_in_vdata;
}

int DataLoader::HardCodedQueueItemOffsetInPdata(int type) {
  if (type == MechanismTypes::
                  kBinReportHelper ||  // &(_nt->_vdata[_ppvar[2*_STRIDE]])
      type == MechanismTypes::kVecStim)
    return 2;
  if (type == MechanismTypes::kALU ||  //&(_nt->_vdata[_ppvar[3*_STRIDE]])
      type == MechanismTypes::kGluSynapse || type == MechanismTypes::kNetStim ||
      type == MechanismTypes::kPatternStim)
    return 3;
  if (type ==
      MechanismTypes::kInhPoissonStim)  //&(_nt->_vdata[_ppvar[6*_STRIDE]])
    return 6;
  return -1;  // means unexistent
}

int DataLoader::HardCodedPPtrOffsetInPdata(int type) {
  if (type == MechanismTypes::kPatternStim || type == MechanismTypes::kALU)
    return 2;
  return -1;  // means unexistent
}

bool DataLoader::HardCodedMechanismHasNoInstances(int type) {
  return (HardCodedMechanismForCoreneuronOnly(type) || type == kVecStim);
}

bool DataLoader::HardCodedMechanismForCoreneuronOnly(int type) {
  // BinReportHelper, BinReport, MemUsage, CoreConfig and ProfileHelper
  // have nodecount 1 but no data and no nodeindices
  return (type == kBinReports || type == kBinReportHelper ||
          type == kCoreConfig || type == kMemUsage || type == kMemUsage ||
          type == kProfileHelper);
}

bool DataLoader::HardCodedEventIsDiscontinuity(Event *e) {
  // This works for the PCP circuit (for now);
  if (e->Type() != EventTypes::kNetCon) return false;

  const int type = ((NetconX *)e)->mech_type_;

  return type == MechanismTypes::kProbAMPANMDA_EMS ||
         type == MechanismTypes::kProbGABAAB_EMS;
}

neuron_id_t DataLoader::GetNeuronIdFromNrnThreadId(int nrn_id) {
  return (neuron_id_t)nrn_threads[nrn_id].presyns[0].gid_;
}

/* execution time for the slowest branch in this locality
 * (for the calculus of max parallelism cap on branch-level parallelism) */
static double slowest_compartment_runtime = 0;

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
  for (NrnThreadMembList *tml = nt->tml; tml != NULL; tml = tml->next) {
    if (HardCodedMechanismHasNoInstances(tml->index)) {
      assert(tml->ml->nodeindices == NULL);
    } else {
      assert(tml->ml->nodeindices != NULL);
    }
    assert(tml->ml->data != NULL);
    assert(nt->ncell == 1);
    data_size_padded += Vectorizer::SizeOf(tml->ml->nodecount) *
                        mechanisms_[mechanisms_map_[tml->index]]->data_size_;
  }

  // map of post- to pre-padding values of pdata
  // only used for branched neurons, otherwise pointers for padded and
  // non-padded layouts are the same
  std::vector<int> data_offsets(
      input_params_->branch_parallelism_ ? data_size_padded : 0, -99999);

  if (input_params_->branch_parallelism_)
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
            || axon_initial_segment_compartments.find(comp->id_) !=
                   axon_initial_segment_compartments.end())  // connection to
                                                             // any AIS
                                                             // compartment
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

  // information about offsets in data and node ids of all instances of all ions
  vector<DataLoader::IonInstancesInfo> ions_instances_info(
      Mechanism::IonTypes::kSizeAllIons);
  int DELETE_TEST = N * 6;
  int DELETE_TEST_nodes = 0;
  // TODO from tqitem we get mech >< instance and which branch it belongs to
  // TODO is tditem of class tqueue.h :: TQItem?
  // Compartment *no_instances_compartment =
  //   new Compartment(-1, -1, -1, -1, -1, -1, -1, -1);
  Compartment *no_instances_compartment = compartments.at(0);
  for (NrnThreadMembList *tml = nt->tml; tml != NULL;
       tml = tml->next)  // For every mechanism
  {
    int type = tml->index;
    Memb_list *ml = tml->ml;  // Mechanisms application to each compartment
    Mechanism *mech = GetMechanismFromType(type);
    assert(mech->type_ == type);
    int ion_offset = mech->GetIonIndex();

    // for every mech instance or compartment this mech is applied to..
    // (Mechs without instances have nodecount 1 and nodeindices NULL)
    DELETE_TEST_nodes += ml->nodecount * mech->data_size_;
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
        // compare before and after offsets
        assert(DELETE_TEST++ == &ml->data[offset_non_padded] - nt->_data);
        // compare nt vs ml offsets
        assert(ml->data[offset_non_padded] ==
               nt->_data[data_total_offset + offset_non_padded]);
        data.push_back(ml->data[offset_non_padded]);
#else
        int offset_padded = Vectorizer::SizeOf(ml->nodecount) * i + n;
        data.push_back(ml->data[offset_padded]);
        assert(ml->data[offset_padded] ==
               nt->_data[data_total_padded_offset + offset_padded]);
        if (input_params_->branch_parallelism_)
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
        int ptype = memb_func[mech->type_].dparam_semantics[i];
        int pd = ml->pdata[pdata_offset_padded];

        // remove extra space added by padding
        //(for pointer to area or ion mech instance)
        if (input_params_->branch_parallelism_ &&
            (ptype == -1 || (ptype > 0 && ptype < 1000))) {
          // if mech has no instances, area is always -1
          if (ptype == -1 /*area*/ && HardCodedMechanismHasNoInstances(type)) {
            assert(pd == -1);
            pdata.push_back(pd);
          } else {
            assert(data_offsets.at(pd) != -99999);
            // offset to non-padded SoA value
            pdata.push_back(data_offsets.at(pd));
          }
        } else
          pdata.push_back(pd);
#endif
      }
      pdata.shrink_to_fit();

      void **vdata = &nt->_vdata[vdata_total_offset];

      bool mech_has_no_instances = HardCodedMechanismHasNoInstances(type);
      /// bool mech_has_no_instances = tml->ml->nodeindices==NULL;
      Compartment *compartment = mech_has_no_instances
                                     ? no_instances_compartment
                                     : compartments.at(ml->nodeindices[n]);
      compartment->AddMechanismInstance(type, n, data.data(), mech->data_size_,
                                        pdata.data(), mech->pdata_size_);
      if (mech->pnt_map_ > 0 || mech->vdata_size_ > 0)  // vdata
      {
        const int point_proc_offset_in_pdata =
            HardCodedPntProcOffsetInPdata(type);
        if (point_proc_offset_in_pdata != -1) {
          const int point_proc_offset_in_vdata =
              HardCodedPntProcOffsetInVdata(type);
          assert(nt->_vdata[pdata[point_proc_offset_in_pdata]] ==
                 &nt->pntprocs[point_proc_total_offset]);
          Point_process *ppn = &nt->pntprocs[point_proc_total_offset];
          Point_process *pp =
              (Point_process *)vdata[point_proc_offset_in_vdata];
          assert(ppn->_i_instance == pp->_i_instance && ppn->_tid == pp->_tid &&
                 ppn->_type == pp->_type);
          compartment->AddSerializedVData((unsigned char *)(void *)pp,
                                          sizeof(Point_process));
          point_proc_total_offset++;
        }

        int rng_offset_in_pdata = HardCodedRNGOffsetInPdata(type);
        if (rng_offset_in_pdata != -1) {
          int rng_offset_in_vdata = HardCodedRNGOffsetInVdata(type);
          assert(!mech_has_no_instances);  // needs to have instances!
          assert(nt->_vdata[pdata[rng_offset_in_pdata]] ==
                 vdata[rng_offset_in_vdata]);
          nrnran123_State *rng = (nrnran123_State *)vdata[rng_offset_in_vdata];
          nrnran123_State *rng2 =
              (nrnran123_State *)nt->_vdata[pdata[rng_offset_in_pdata]];
          assert(rng->c == rng2->c && rng->r == rng2->r &&
                 rng->which_ == rng2->which_);
          assert(rng != NULL);
          compartment->AddSerializedVData((unsigned char *)(void *)rng,
                                          sizeof(nrnran123_State));
        }

        int queue_item_offset_in_pdata = HardCodedQueueItemOffsetInPdata(type);
        if (queue_item_offset_in_pdata != -1) {
          // TODO copy value not pointer
          void *q_item_ptr = nt->_vdata[pdata[queue_item_offset_in_pdata]];
          if (q_item_ptr == nullptr) {
            if (!silence_queue_item_warning) {
              fprintf(stderr,
                      "Warning: Locality %d, mech %d has NULL TQItem (void*)\n",
                      wrappers::MyRankId(), type);
              silence_queue_item_warning = true;
            }
            // work around, put something which is 8 bytes long
            q_item_ptr = new unsigned char[sizeof(void *)];
            memset(q_item_ptr, 0, sizeof(void *));
          }
          compartment->AddSerializedVData((unsigned char *)(void *)q_item_ptr,
                                          sizeof(void *));
        }

        int p_ptr_offset_in_pdata = HardCodedPPtrOffsetInPdata(type);
        if (p_ptr_offset_in_pdata != -1) {
          // TODO copy value not pointer
          void *p_ptr = nt->_vdata[pdata[p_ptr_offset_in_pdata]];
          if (p_ptr == nullptr) {
            if (!silence_p_ptr_warning) {
              fprintf(stderr,
                      "Warning: Locality %d, mech %d has NULL p_ptr (void*)\n",
                      wrappers::MyRankId(), type);
              silence_p_ptr_warning = true;
            }
            // work around, put something which is 8 bytes long
            p_ptr = new unsigned char[sizeof(void *)];
            memset(p_ptr, 0, sizeof(void *));
          }
          compartment->AddSerializedVData((unsigned char *)(void *)p_ptr,
                                          sizeof(void *));
        }

        vdata_total_offset += (unsigned)mech->vdata_size_;
      }
    }
    data_total_offset += mech->data_size_ * ml->nodecount;
    data_total_padded_offset +=
        mech->data_size_ * Vectorizer::SizeOf(ml->nodecount);
  }
  assert(nt->_nvdata == vdata_total_offset);

  // assert(DELETE_TEST == nt->_ndata);
  // assert(DELETE_TEST_nodes == nt->_ndata - N * 6);
  for (Compartment *comp : compartments) comp->ShrinkToFit();

  data_offsets.clear();

  //======= 3 - reconstruct NetCons =====================
  map<neuron_id_t, vector<NetconX *>> netcons;  // netcons per pre-synaptic
                                                // neuron id)
  for (int n = 0; n < nt->n_netcon; ++n) {
    NetCon *nc = nt->netcons + n;
    assert(netcon_srcgid.size() >
           0);  // if size==0, then setup_cleanup() in nrn_setup.cpp was called
    assert(nt->id < netcon_srcgid.size());
    assert(netcon_srcgid.at(nt->id) != NULL);
    int srcgid = netcon_srcgid[nt->id][n];
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

  // add mechs instances without compartment to end of compartments (if not
  // soma)
  if (no_instances_compartment != compartments.at(0))
    compartments.push_back(no_instances_compartment);
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
  std::vector<ReportConfiguration> configs;
  nrn_init_and_load_data(argc, argv, &configs, nrnmpi_under_nrncontrol,
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

  // TODO has to ble cleared at some point
  Mechanism::time_spent_in_mechs_mutex_ = hpx_lco_sema_new(1);

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

    if (input_params_->graph_mechs_parallelism_) {
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
  if (neurox::ParallelExecution() && input_params_->branch_parallelism_) {
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
  NEUROX_MEM_UNPIN;
}

hpx_action_t DataLoader::Init = 0;
int DataLoader::Init_handler() {
  NEUROX_MEM_PIN(uint64_t);
  all_neurons_gids_ = new std::vector<int>();
  locality_mutex_ = hpx_lco_sema_new(1);

  Statistics::CommCount::counts.point_to_point_count = 0;
  Statistics::CommCount::counts.reduce_count = 0;
  Statistics::CommCount::counts.spike_count = 0;

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
  if (input_params_->locality_comm_reduce_ ||
      synchronizer_->LocalitySyncInterval() == -1) {
    assert(locality::netcons_branches_ == nullptr);
    locality::neurons_ = new vector<hpx_t>();
    locality::netcons_branches_ = new map<neuron_id_t, vector<hpx_t>>();
    locality::netcons_somas_ = new map<neuron_id_t, vector<hpx_t>>();
#if defined(PRINT_TIME_DEPENDENCY) or defined(PRINT_TIME_DEPENDENCY_MUTEX) or \
    defined(PRINT_TIME_DEPENDENCY_STEP_SIZE)
    locality::from_hpx_to_gid = new map<hpx_t, neuron_id_t>();
#endif
  }

  NEUROX_MEM_UNPIN;
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
    if (input_params_->locality_comm_reduce_ || input_params_->scheduler_) {
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
      if (input_params_->graph_mechs_parallelism_) {
        Mechanism *mech = neurox::mechanisms_[m];
        for (int d = 0; d < mech->dependencies_count_; d++) {
          Mechanism *parent = GetMechanismFromType(mech->dependencies_[d]);
          if (strcmp("SK_E2", mech->memb_func_.sym) == 0 &&
              strcmp("ca_ion", parent->memb_func_.sym) == 0)
            continue;  // TODO hard-coded exception

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
    fprintf(file_mechs, "digraph G\n{ bgcolor=%s; %s\n", "transparent",
            !input_params_->graph_mechs_parallelism_
                ? "layout=circo; scale=0.23;"
                : "");
    fprintf(file_mechs, "graph [ratio=0.3];\n");
    fprintf(file_mechs, "%s [style=filled, shape=Mdiamond, fillcolor=beige];\n",
            "start");
    fprintf(file_mechs, "%s [style=filled, shape=Mdiamond, fillcolor=beige];\n",
            "end");
    fprintf(file_mechs, "\"%s (%d)\" [style=filled, fillcolor=beige];\n",
            GetMechanismFromType(CAP)->memb_func_.sym, CAP);
    if (!input_params_->graph_mechs_parallelism_) {
      fprintf(file_mechs, "end -> start [color=transparent];\n");
      fprintf(file_mechs, "start -> \"%s (%d)\";\n",
              GetMechanismFromType(CAP)->memb_func_.sym, CAP);
    }
    for (int m = 0; m < mechanisms_count_; m++) {
      Mechanism *mech = neurox::mechanisms_[m];

      if (mech->pnt_map_ > 0)  // if is point process make it dotted
        fprintf(file_mechs, "\"%s (%d)\" [style=dashed%s];\n",
                mech->memb_func_.sym, mech->type_,
                HardCodedMechanismForCoreneuronOnly(mech->type_)
                    ? ", color=gray, fontcolor=gray"
                    : "");

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

      if (input_params_->graph_mechs_parallelism_)
        for (int d = 0; d < mech->dependencies_count_; d++) {
          Mechanism *parent = GetMechanismFromType(mech->dependencies_[d]);
          if (strcmp("SK_E2", mech->memb_func_.sym) == 0 &&
              strcmp("ca_ion", parent->memb_func_.sym) == 0)
            continue;  // TODO: hard-coded exception
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
          "vars, noInstances %d\n",
          mech->memb_func_.sym, mech->type_, mech->data_size_,
          mech->pdata_size_, mech->is_artificial_, mech->pnt_map_,
          mech->is_ion_, mech->sym_length_, mech->successors_count_,
          mech->dependencies_count_,
          mech->state_vars_ ? mech->state_vars_->count_ : 0,
          HardCodedMechanismHasNoInstances(mech->type_));
    }
  }
#endif

  all_neurons_gids_->clear();
  delete all_neurons_gids_;
  all_neurons_gids_ = nullptr;

  nrn_setup_cleanup();

#if defined(NDEBUG)
  // if not on debug, there's no CoreNeuron comparison, so data can be
  // cleaned-up now, except global_ion_map
  DataLoader::CleanCoreneuronData(false);
#endif

  if (input_params_->output_statistics_)
    tools::LoadBalancing::LoadBalancing::PrintLoadBalancingTable();

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

  hpx_lco_delete_sync(locality_mutex_);
  return neurox::wrappers::MemoryUnpin(target);
}

void DataLoader::GetSubSectionFromCompartment(deque<Compartment *> &sub_section,
                                              Compartment *top_comp) {
  // parent added first, to respect solvers logic, of parents' ids first
  sub_section.push_back(top_comp);
  for (int c = 0; c < top_comp->branches_.size(); c++)
    GetSubSectionFromCompartment(sub_section, top_comp->branches_.at(c));
}

void DataLoader::GetMechInstanceMap(const deque<Compartment *> &compartments,
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
      ppi.mech_instance = (offset_t)(
          *mech_instances_map)[neurox::mechanisms_map_[ppi.mech_type]]
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
  for (const auto comp : compartments)
    if (comp->id_ != -1) data.push_back(comp->rhs_);
  for (const auto comp : compartments)
    if (comp->id_ != -1) data.push_back(comp->d_);
  for (const auto comp : compartments)
    if (comp->id_ != -1) data.push_back(comp->a_);
  for (const auto comp : compartments)
    if (comp->id_ != -1) data.push_back(comp->b_);
  for (const auto comp : compartments)
    if (comp->id_ != -1) data.push_back(comp->v_);
  for (const auto comp : compartments)
    if (comp->id_ != -1) data.push_back(comp->area_);
  for (const auto comp : compartments)
    if (comp->id_ != -1) p.push_back(comp->p_);
  assert(p.size() > 0);  // zero-sized sub-section

  ////// Tree of neurons: convert from neuron- to branch-level
  std::map<int, int> from_old_to_new_compartment_id;
  for (const Compartment *comp : compartments)
    if (comp->id_ != -1) from_old_to_new_compartment_id[comp->id_] = n++;

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
      int new_id =
          comp->id_ == -1 ? -1 : from_old_to_new_compartment_id.at(comp->id_);
      nodes_indices_mechs[mech_offset].push_back(new_id);

      instances_count[mech_offset]++;
      comp_data_offset += mech->data_size_;
      comp_pdata_offset += mech->pdata_size_;

      if (mech->pnt_map_ > 0 || mech->vdata_size_ > 0) {
        size_t total_vdata_size = HardCodedVdataSize(mech->type_);
        for (int vi = comp_vdata_offset;
             vi < comp_vdata_offset + total_vdata_size; vi++) {
          assert(vi < comp->vdata_serialized_.size());
          vdata_mechs[mech_offset].push_back(comp->vdata_serialized_.at(vi));
        }
        /* //TODO put it back
        vdata_mechs[mech_offset].insert(
            vdata_mechs[mech_offset].end(), &comp->vdata_[comp_vdata_offset],
            &comp->vdata_[comp_vdata_offset + total_vdata_size]);
        */
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
    int mech_has_no_instances = HardCodedMechanismHasNoInstances(mech->type_);

    // for all instances
    for (int i = 0; i < nodes_indices_mechs[m].size(); i++) {
      if (!HardCodedMechanismHasNoInstances(mech->type_))
        nodes_indices.push_back(nodes_indices_mechs[m][i]);

      // pdata calculation for branch-parallelism
      assert(ion_instance_to_data_offset.find(
                 make_pair(mech->type_, nodes_indices_mechs[m][i])) ==
             ion_instance_to_data_offset.end());

      // set look-up map with beginning of data instance offset
      if (mech->is_ion_ && input_params_->branch_parallelism_) {
        assert(!mech_has_no_instances);
        ion_instance_to_data_offset[make_pair(
            mech->type_, nodes_indices_mechs[m][i])] = data.size();
      }

      // insert data (TODO this was before the previous condition)
      data.insert(data.end(), &data_mechs[m][data_offset],
                  &data_mechs[m][data_offset + mech->data_size_]);
      // TODO de we need to account from SIMD spacing here???
      data_offset += mech->data_size_;

      if (mech->pnt_map_ > 0 || mech->vdata_size_ > 0) {
        int total_vdata_size = HardCodedVdataSize(mech->type_);
        assert(total_vdata_size > 0);
        vdata.insert(vdata.end(), &vdata_mechs[m][vdata_offset],
                     &vdata_mechs[m][vdata_offset + total_vdata_size]);
        vdata_offset += total_vdata_size;
      }

      // if we need to recalculate offsets or remove padding
      if (input_params_->branch_parallelism_) {
        for (int p = pdata_offset; p < pdata_offset + mech->pdata_size_; p++) {
          offset_t pd = pdata_mechs.at(m).at(p);
          int ptype = memb_func[mech->type_].dparam_semantics[p - pdata_offset];
          switch (ptype) {
            case -1:  //"area" (6th field)
            {
              if (HardCodedMechanismHasNoInstances(mech->type_)) {
                assert(pd == -1);  // no instances, no compartment, no area!
                pdata_mechs.at(m).at(p) = pd;
              } else {
                assert(pd >= N * 5 && pd < N * 6);
                offset_t old_id = pd - N * 5;
                offset_t new_id = from_old_to_new_compartment_id.at(old_id);
                pdata_mechs.at(m).at(p) = n * 5 + new_id;
              }
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
                // assert(node_id == new_node_id); //only for branch-parallelism
                // with single section
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

double DataLoader::BenchmarkSubSection(
    int N, const deque<Compartment *> &subsection,
    vector<DataLoader::IonInstancesInfo> &ions_instances_info) {
  // common vars to all compartments
  int dumb_threshold_offset = 0;
  int nrn_threadId = -1;
  hpx_t soma_branch_addr = HPX_NULL;

  // timing vars
  double time_elapsed = 0;

  // serialization of neurons
  offset_t n;              // number of compartments in branch
  vector<floble_t> data;   // compartments info (RHS, D, A, B, V, AREA)*n
  vector<offset_t> pdata;  // pointers to data
  vector<offset_t> p;      // parent nodes index
  vector<offset_t> instances_count(mechanisms_count_);
  vector<offset_t> nodes_indices;
  vector<unsigned char> vdata;  // Serialized Point Processes and Random123
  vector<floble_t> vecplay_t;
  vector<floble_t> vecplay_y;
  vector<PointProcInfo> vecplay_info;
  vector<NetconX> branch_netcons;
  vector<neuron_id_t> branch_netcons_pre_id;
  vector<floble_t> branch_weights;

  // allocate GAS memory for this temporary branch
  hpx_t temp_branch_addr =
      hpx_gas_alloc_local(1, sizeof(Branch), Vectorizer::kMemoryAlignment);

  // mech-offset -> ( map[old instance]->to new instance )
  vector<map<int, int>> mech_instances_map(neurox::mechanisms_count_);
  GetMechInstanceMap(subsection, mech_instances_map);

  n = GetBranchData(subsection, data, pdata, vdata, p, instances_count,
                    nodes_indices, N, ions_instances_info, &mech_instances_map);
  GetVecPlayBranchData(subsection, vecplay_t, vecplay_y, vecplay_info,
                       &mech_instances_map);
  GetNetConsBranchData(subsection, branch_netcons, branch_netcons_pre_id,
                       branch_weights, &mech_instances_map);

  // initialize datatypes and graph-parallelism shadow vecs offsets
  neuron_id_t soma_gid = 0;
  floble_t soma_ap_threshold = 999;

  hpx_call_sync(
      temp_branch_addr, Branch::Init, &time_elapsed,
      sizeof(time_elapsed),  // output
      &n, sizeof(offset_t), &nrn_threadId, sizeof(int), &dumb_threshold_offset,
      sizeof(int), data.size() > 0 ? data.data() : nullptr,
      sizeof(floble_t) * data.size(), pdata.size() > 0 ? pdata.data() : nullptr,
      sizeof(offset_t) * pdata.size(), instances_count.data(),
      instances_count.size() * sizeof(offset_t), nodes_indices.data(),
      nodes_indices.size() * sizeof(offset_t), &soma_branch_addr, sizeof(hpx_t),
      nullptr, 0,                             // no branches
      p.data(), sizeof(offset_t) * p.size(),  // force use of parent index
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
      sizeof(unsigned char) * vdata.size(), &soma_gid, sizeof(neuron_id_t),
      &soma_ap_threshold, sizeof(floble_t));

  /* this mem pin works because benchmark neurons are allocated locally*/
  Branch *branch = NULL;
  int err = hpx_gas_try_pin(temp_branch_addr, (void **)&branch);
  assert(err != 0);

  // initialize datatypes and graph-parallelism shadow vecs offsets
  interpolators::BackwardEuler::Finitialize2(branch);

  // benchmark execution time of 10 times 0.1msec
  time_elapsed = 0;
  hpx_time_t now = hpx_time_now();
  if (input_params_->interpolator_ ==
      interpolators::InterpolatorIds::kBackwardEuler) {
    const int comm_steps =
        (neurox::min_synaptic_delay_ + 0.00001) / input_params_->dt_;
    for (int i = 0; i < 20; i++) {
      for (int i = 0; i < comm_steps; i++)
        interpolators::BackwardEuler::Step(branch, true);
    }
  } else  // variable time step
  {
    for (int i = 1; i <= 20; i++) {
      // TODO not implemented
      // interpolators::VariableTimeStep::StepTo(branch,0.1*i);
    }
  }
  time_elapsed += hpx_time_elapsed_ms(now) / 1e3;
  time_elapsed /= 20;

  hpx_call_sync(temp_branch_addr, Branch::Clear, nullptr, 0);
  hpx_gas_unpin(temp_branch_addr);
  hpx_gas_clear_affinity(temp_branch_addr);
  return time_elapsed;
}

int DataLoader::GetNumberOfInstanceCompartments(
    const deque<Compartment *> &compartments) {
  for (Compartment *comp : compartments)
    if (comp->id_ == -1) return compartments.size() - 1;
  return compartments.size();
}

hpx_t DataLoader::CreateBranch(
    const int nrn_threadId, hpx_t soma_branch_addr,
    const deque<Compartment *> &all_compartments, Compartment *top_compartment,
    vector<DataLoader::IonInstancesInfo> &ions_instances_info,
    double neuron_runtime, int thvar_index /*AIS*/,
    floble_t ap_threshold /*AIS*/, int assigned_locality) {
  // remove no_instances_mechanism from end of all_compartments
  int N = GetNumberOfInstanceCompartments(all_compartments);
  bool is_soma = soma_branch_addr == HPX_NULL;
  bool is_AIS = thvar_index != -1;

  assert(top_compartment != NULL);
  offset_t n;              // number of compartments in branch
  vector<floble_t> data;   // compartments info (RHS, D, A, B, V, AREA)*n
  vector<offset_t> pdata;  // pointers to data
  vector<offset_t> p;      // parent nodes index
  vector<offset_t> instances_count(mechanisms_count_);
  vector<offset_t> nodes_indices;
  vector<hpx_t> branches;
  vector<unsigned char> vdata;  // Serialized Point Processes and Random123

  // Vector Play instances
  vector<floble_t> vecplay_t;
  vector<floble_t> vecplay_y;
  vector<PointProcInfo> vecplay_info;

  // branch NetCons
  vector<NetconX> branch_netcons;
  vector<neuron_id_t> branch_netcons_pre_id;
  vector<floble_t> branch_weights;

  Compartment *bottom_compartment = nullptr;
  deque<Compartment *> subsection;
  double subsection_runtime = -1;

  /* defult value (-1) means use this locality */
  if (assigned_locality < 0) assigned_locality = hpx_get_my_rank();

  /* get the runtime of SIMD structured neuron */
  if (is_soma /*first run*/) {
    // Exclude no instances mechanisms from neuron time
    // TODO this is the right version, after fixing TODO below
    // GetSubSectionFromCompartment(subsection, all_compartments.at(0));
    /* compute total runtime of neuron */
    // neuron_runtime =
    //       BenchmarkSubSection(N, subsection, ions_instances_info, true,
    //       false);

    /* compute total runtime of neuron */
    neuron_runtime =
        input_params_->interpolator_ !=
                interpolators::InterpolatorIds::kBackwardEuler
            ? -1
            :  // disable for CVODE
            BenchmarkSubSection(N, all_compartments, ions_instances_info);
  }

  /* No branch-parallelism (a la Coreneuron) */
  if (!input_params_->branch_parallelism_) {
    n = GetBranchData(all_compartments, data, pdata, vdata, p, instances_count,
                      nodes_indices, N, ions_instances_info, NULL);
    GetVecPlayBranchData(all_compartments, vecplay_t, vecplay_y, vecplay_info,
                         NULL);
    GetNetConsBranchData(all_compartments, branch_netcons,
                         branch_netcons_pre_id, branch_weights, NULL);

    /* For querying the table we use neuron's exec time */
    subsection = all_compartments;
    subsection_runtime = neuron_runtime;
    if (input_params_->load_balancing_) {
      /* assign branch locally*/
      hpx_call_sync(HPX_THERE(0), tools::LoadBalancing::QueryLoadBalancingTable,
                    &assigned_locality, sizeof(int));  // output
    }
    /* branch - parallelism */
  } else {
    /* estimated max execution time per subregion for this neuron*/
    double max_work_per_subtree = LoadBalancing::GetMaxWorkPerBranchSubTree(
        neuron_runtime, GetMyNrnThreadsCount());

#ifndef NDEBUG
    if (is_soma)
      printf(
          "=== nrn_id %d: neuron time %.5f ms, max runtime per subtree %.5f "
          "ms\n",
          nrn_threadId, neuron_runtime, max_work_per_subtree);
#endif

    /*pack compartments until finding bifurcation */
    /* new subsection: iterate until bifurcation or max time is reached */

    // if remaining arborization fits, use it as a subsection
    GetSubSectionFromCompartment(subsection, top_compartment);

    // This is just necessary to debug branch parallelism with single sections
    // against no branch-parallelism
    std::sort(subsection.begin(), subsection.end(), CompareCompartmentPtrId);

    subsection_runtime =
        BenchmarkSubSection(N, subsection, ions_instances_info);
    // printf("TEST: arborization starting on comp %d, time %f\n",
    // top_compartment->id_, subsection_runtime);

    // if not: pack compartments until it does, or find a bifurcation
    if (subsection_runtime > max_work_per_subtree) {
      subsection.clear();
      subsection.push_back(top_compartment);
      subsection_runtime = -1;
      for (bottom_compartment = top_compartment;
           bottom_compartment->branches_.size() == 1;
           bottom_compartment = bottom_compartment->branches_.front()) {
        // test current subsection
        subsection_runtime =
            BenchmarkSubSection(N, subsection, ions_instances_info);
        // printf("TEST: Subsections starting on comp %d, length %d, time %f vs
        // %f\n",
        //       bottom_compartment->id_, subsection.size(), subsection_runtime,
        //       max_work_per_subtree);

        // keeps adding compartments to subsection until remaining tree fits in
        // runtime
        if (subsection_runtime < max_work_per_subtree)
          // add a single compartment to subsection and re-try
          subsection.push_back(bottom_compartment->branches_.front());
        else
          break;
      }
      // three conditions that neets to be verified for this to work
      assert(subsection_runtime != -1 && subsection.size() > 0 &&
             bottom_compartment);
    }

    //#ifndef NDEBUG
    printf("--- nrn_id %d: %s %d, length %d, runtime %.5f ms\n", nrn_threadId,
           is_soma ? "soma" : (is_AIS ? "AIS" : "subsection"),
           top_compartment->id_, subsection.size(), subsection_runtime);
    //#endif

    /* create serialized sub-section from compartments in subsection*/
    // mech-offset -> ( map[old instance]->to new instance )
    vector<map<int, int>> mech_instances_map(neurox::mechanisms_count_);
    GetMechInstanceMap(subsection, mech_instances_map);

    n = GetBranchData(subsection, data, pdata, vdata, p, instances_count,
                      nodes_indices, N, ions_instances_info,
                      &mech_instances_map);
    GetVecPlayBranchData(subsection, vecplay_t, vecplay_y, vecplay_info,
                         &mech_instances_map);
    GetNetConsBranchData(subsection, branch_netcons, branch_netcons_pre_id,
                         branch_weights, &mech_instances_map);

    // Load balancing matters here:
    if (input_params_->load_balancing_) {
      /* estimated max executiom time allowed per locality */
      double max_work_per_subsection =
          LoadBalancing::GetMaxWorkPerBranchSubSection(neuron_runtime,
                                                       GetMyNrnThreadsCount());

      assert(0);  // TODO double check
                  /* TODO: check: soma and AIS cant be split due to AP threshold
                   * communication */

      //#ifndef NDEBUG
      if (is_soma)
        printf("==== neuron time %.5f ms, max runtime per locality %.5fms\n",
               neuron_runtime, max_work_per_subsection);
      //#endif

      /* if subsection_complexity == 0, keep all subtrees in the
       * same locality (get rank of soma, or use previously assigned
       * assigned_locality otherwise) */
      if ((max_work_per_subsection == 0 && is_soma)
          /* assign remaining arborization to different locality if it fits in
           * the max runtime per locality and is not too small (to reduce number
           * of remote small branches)*/
          || (max_work_per_subsection > 0 &&
              assigned_locality == hpx_get_my_rank() &&
              subsection_runtime < max_work_per_subsection)) {
        // subsection_runtime < max_work_per_locality &&
        // subsection_runtime > max_work_per_locality * 0.5) {
        /* ask master rank where to allocate this arborization, update table*/
        hpx_call_sync(HPX_THERE(0),
                      tools::LoadBalancing::QueryLoadBalancingTable,  /*action*/
                      &assigned_locality, sizeof(assigned_locality)); /*output*/
      }
    }
  }

  if (input_params_->output_statistics_) {
    /* for the LPT table and statistics, we use the final data struct runtime*/
    if (subsection_runtime == -1)
      subsection_runtime =
          BenchmarkSubSection(N, subsection, ions_instances_info);

    /* tell master rank to update entry in Least-Processing-Time table */
    hpx_call_sync(HPX_THERE(0), tools::LoadBalancing::UpdateLoadBalancingTable,
                  nullptr, 0,                           // output
                  &subsection_runtime, sizeof(double),  // input[0]
                  &assigned_locality, sizeof(int));     // input[1]

    //#ifndef NDEBUG
    printf(
        "- %s %d, length %d, nrn_id %d, actual runtime %.6f ms, allocated to "
        "%s rank "
        "%d\n",
        is_soma ? "soma" : (is_AIS ? "AIS" : "subsection"),
        top_compartment->id_, n, nrn_threadId, subsection_runtime,
        assigned_locality == hpx_get_my_rank() ? "local" : "remote",
        assigned_locality);
    //#endif
  }

  /* allocate GAS mem for subsection in the appropriate locality */
  hpx_t branch_addr = hpx_gas_alloc_local_at_sync(1, sizeof(Branch),
                                                  Vectorizer::kMemoryAlignment,
                                                  HPX_THERE(assigned_locality));

  // update hpx address of soma
  soma_branch_addr = is_soma ? branch_addr : soma_branch_addr;

  /* allocate children branches recursively (if any) */
  if (bottom_compartment) {
    // TODO   make sure soma is not split! or this AIS test fails!!!!
    for (size_t c = 0; c < bottom_compartment->branches_.size(); c++)
      branches.push_back(CreateBranch(
          nrn_threadId, soma_branch_addr, all_compartments,
          bottom_compartment->branches_[c], ions_instances_info, neuron_runtime,
          is_soma && c == 0 ? thvar_index - n : -1, assigned_locality));
    /*offset in AIS = offset in soma - nt->end */

    if (is_soma) thvar_index = -1;
  }

  /* initialize subsection on the appropriate locality */
  neuron_id_t neuron_id =
      is_soma ? GetNeuronIdFromNrnThreadId(nrn_threadId) : -1;
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
      sizeof(unsigned char) * vdata.size(), &neuron_id, sizeof(neuron_id_t),
      &ap_threshold, sizeof(floble_t));

  if (is_soma) {
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

  // reconstruct map of locality to branch netcons (if needed)
  hpx_t top_branch_addr =
      local->soma_ ? target : local->branch_tree_->top_branch_addr_;
  if (input_params_->locality_comm_reduce_) {
    hpx_lco_sema_p(input::DataLoader::locality_mutex_);
    for (auto &nc : local->netcons_) {
      const neuron_id_t pre_neuron_id = nc.first;
      (*locality::netcons_branches_)[pre_neuron_id].push_back(target);
      (*locality::netcons_somas_)[pre_neuron_id].push_back(top_branch_addr);
      // duplicates will be deleted in DataLoader::Finalize
    }
    hpx_lco_sema_v_sync(input::DataLoader::locality_mutex_);
  }

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
        if (input_params_->synchronizer_ == SynchronizerIds::kTimeDependency) {
          dependencies.push_back(make_pair(src_gid, min_delay));
        }
      }
    }
  }

  // inform pre-syn neuron that he connects to this branch or locality
  // (repeated entries will be removed by DataLoader::FilterLocalitySynapses)
  hpx_t netcons_lco = hpx_lco_and_new(netcons.size());
  int my_gid = local->soma_ ? local->soma_->gid_ : -1;
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

hpx_action_t DataLoader::FilterRepeatedAndLinearizeContainers = 0;
int DataLoader::FilterRepeatedAndLinearizeContainers_handler() {
  NEUROX_MEM_PIN(Branch);

  if (input_params_->locality_comm_reduce_) {
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
  }

  // convert synapses to linear synapses representation
  if (input_params_->linearize_containers_) local->soma_->LinearizeContainers();

  NEUROX_MEM_UNPIN
}

void DataLoader::RegisterHpxActions() {
  wrappers::RegisterZeroVarAction(
      DataLoader::FilterRepeatedAndLinearizeContainers,
      DataLoader::FilterRepeatedAndLinearizeContainers_handler);
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
