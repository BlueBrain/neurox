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
#include "coreneuron/nrnoc/nrnoc_decl.h"
#include "coreneuron/nrnoc/nrnoc_decl.h"  //nrn_is_ion()
#include "coreneuron/utils/memory_utils.h"
#include "coreneuron/utils/randoms/nrnran123.h"  //RNG data structures

using namespace std;
using namespace neurox::input;
using namespace neurox::algorithms;

FILE *DataLoader::fileNetcons = nullptr;
hpx_t DataLoader::all_neurons_mutex = HPX_NULL;
std::vector<hpx_t> *DataLoader::my_neurons_addr = nullptr;
std::vector<int> *DataLoader::my_neurons_gids = nullptr;
std::vector<int> *DataLoader::all_neurons_gids = nullptr;
tools::LoadBalancing *DataLoader::loadBalancing = nullptr;

neuron_id_t DataLoader::GetNeuronIdFromNrnThreadId(int nrn_id) {
  return (neuron_id_t)nrn_threads[nrn_id].presyns[0].gid_;
}

void DataLoader::PrintSubClustersToFile(FILE *fileCompartments,
                                        Compartment *topCompartment) {
  if (input_params->outputCompartmentsDot) {
    assert(topCompartment != NULL);
    fprintf(fileCompartments, "subgraph cluster_%d { ", topCompartment->id);
    if (topCompartment->id == 0)
      fprintf(fileCompartments, "label=\"SOMA\"; ");
    else if (topCompartment->id == 2)
      fprintf(fileCompartments, "label=\"AIS\"; ");
    Compartment *comp = NULL;
    for (comp = topCompartment; comp->branches.size() == 1;
         comp = comp->branches.front())
      fprintf(fileCompartments, "%d; ", comp->id);
    fprintf(fileCompartments, "%d };\n", comp->id);
    for (int c = 0; c < comp->branches.size(); c++)
      PrintSubClustersToFile(fileCompartments, comp->branches[c]);
  }
}

PointProcInfo DataLoader::GetPointProcInfoFromDataPointer(NrnThread *nt,
                                                          double *pd,
                                                          size_t size) {
  PointProcInfo ppi;
  ppi.nodeId = -1;
  ppi.size = size;
  bool found = false;
  for (NrnThreadMembList *tml = nt->tml; !found && tml != NULL;
       tml = tml->next)  // For every mechanism
  {
    int type = tml->index;
    Mechanism *mech = GetMechanismFromType(type);
    Memb_list *ml = tml->ml;
    for (int n = 0; n < ml->nodecount; n++)
      for (int i = 0; i < mech->dataSize; i++) {
#if LAYOUT == 1
        int dataOffset = mech->dataSize * n + i;
#else
        int dataOffset = tools::Vectorizer::SizeOf(ml->nodecount) * i + n;
#endif
        if (&ml->data[dataOffset] != pd) continue;  // if not this variable

        ppi.nodeId = ml->nodeindices[n];
        ppi.mechType = type;
        ppi.mechInstance = (offset_t)n;
        ppi.instanceDataOffset = i;
        found = true;
        break;
      }
  }
  assert(found);
  return ppi;
}

int DataLoader::CreateNeuron(int neuron_idx, void *) {
  NrnThread *nt = &nrn_threads[neuron_idx];
  neuron_id_t neuronId = GetNeuronIdFromNrnThreadId(nt->id);
  int N = nt->end;

  // if data is permuted, this method fails.
  assert(!use_interleave_permute && !use_solve_interleave &&
         nt->_permute == NULL);

  // map of padded to non-padded offsets of data
  size_t dataSizePadded = 6 * tools::Vectorizer::SizeOf(N);
  for (NrnThreadMembList *tml = nt->tml; tml != NULL; tml = tml->next)
    dataSizePadded += tools::Vectorizer::SizeOf(tml->ml->nodecount) *
                      mechanisms[mechanisms_map[tml->index]]->dataSize;

  // map of post- to pre-padding values of pdata
  // only used for branched neurons, otherwise pointers for padded and
  // non-padded layouts are the same
  std::vector<int> dataOffsets(
      input_params->branchingDepth > 0 ? dataSizePadded : 0, -99999);

  if (input_params->branchingDepth > 0)
    for (int n = 0; n < N; n++)
      for (int i = 0; i < 6; i++) {
        int offsetPadded = tools::Vectorizer::SizeOf(N) * i + n;
        int offsetNonPadded = N * i + n;
        dataOffsets[offsetPadded] = offsetNonPadded;
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
    Compartment *parentCompartment = compartments.at(nt->_v_parent_index[n]);
    parentCompartment->AddChild(compartments.at(n));
  }

  if (input_params->outputCompartmentsDot) {
    set<int> axonInitSegmentCompartments;
    FILE *fileCompartments = fopen(
        string("compartments_" + to_string(neuronId) + ".dot").c_str(), "wt");
    fprintf(fileCompartments, "graph G_%d\n{ bgcolor=%s; \n", neuronId,
            "transparent");
    fprintf(fileCompartments,
            "graph [fontname=helvetica, style=filled, color=blue, "
            "fillcolor=floralwhite];\n");
    fprintf(fileCompartments,
            "node [fontname=helvetica, shape=cylinder, color=gray, "
            "style=filled, fillcolor=white];\n");
    fprintf(fileCompartments, "edge [fontname=helvetica, color=gray];\n");
    PrintSubClustersToFile(fileCompartments,
                           compartments.at(0));  // add subclusters
    for (Compartment *comp : compartments)       // draw edges
      for (int c = 0; c < comp->branches.size(); c++) {
        bool isSoma = comp->id == 1;  // bottom of soma
        Compartment *child = comp->branches.at(c);
        if ((isSoma && c == 0)  // connection from some to AIS
            ||
            axonInitSegmentCompartments.find(comp->id) !=
                axonInitSegmentCompartments
                    .end())  // connection to any AIS compartment
        {
          // reverse connection, so that it plots AIS on top of soma in dot file
          fprintf(fileCompartments, "%d -- %d%s;\n", child->id, comp->id,
                  comp->branches.size() == 1 ? "" : " [color=black]");
          axonInitSegmentCompartments.insert(child->id);
        } else
          fprintf(fileCompartments, "%d -- %d%s;\n", comp->id, child->id,
                  comp->branches.size() == 1 ? "" : " [color=black]");
      }
    fprintf(fileCompartments, "}\n");
    fclose(fileCompartments);
  }

  //======= 2 - reconstructs mechanisms instances ========
  unsigned vdataTotalOffset = 0;
  unsigned dataTotalOffset = N * 6;  // no padding
  unsigned dataTotalPaddedOffset =
      tools::Vectorizer::SizeOf(N) * 6;  // with padding
  unsigned pointProcTotalOffset = 0;

  // information about offsets in data and node ifs of all instances of all ions
  vector<DataLoader::IonInstancesInfo> ionsInstancesInfo(
      Mechanism::IonTypes::kSizeAllIons);
  for (NrnThreadMembList *tml = nt->tml; tml != NULL;
       tml = tml->next)  // For every mechanism
  {
    int type = tml->index;
    Memb_list *ml = tml->ml;  // Mechanisms application to each compartment
    Mechanism *mech = GetMechanismFromType(type);
    assert(mech->type == type);
    int ionOffset = mech->GetIonIndex();
    for (int n = 0; n < ml->nodecount; n++)  // for every mech instance (or
                                             // compartment this mech is applied
                                             // to)
    {
      if (mech->isIon) {
        if (n == 0) {
          ionsInstancesInfo[ionOffset].mechType = type;
          ionsInstancesInfo[ionOffset].dataStart = dataTotalOffset;
          ionsInstancesInfo[ionOffset].dataEnd =
              dataTotalOffset + ml->nodecount * mech->dataSize;
        }
        ionsInstancesInfo[ionOffset].nodeIds.push_back(ml->nodeindices[n]);
      }

      std::vector<double> data;
      for (int i = 0; i < mech->dataSize; i++) {
        int offsetNonPadded = mech->dataSize * n + i;
#if LAYOUT == 1
        data.push_back(ml->data[offsetNonPadded]);
        assert(ml->data[offsetNonPadded] ==
               nt->_data[dataTotalOffset + offsetNonPadded]);
#else
        int offsetPadded = tools::Vectorizer::SizeOf(ml->nodecount) * i + n;
        data.push_back(ml->data[offsetPadded]);
        assert(ml->data[offsetPadded] ==
               nt->_data[dataTotalPaddedOffset + offsetPadded]);
        if (input_params->branchingDepth > 0)
          dataOffsets[dataTotalPaddedOffset + offsetPadded] =
              dataTotalOffset + offsetNonPadded;
#endif
      }
      data.shrink_to_fit();

      std::vector<int> pdata;
      for (int i = 0; i < mech->pdataSize; i++) {
#if LAYOUT == 1
        int pdataOffsetNonPadded = mech->pdataSize * n + i;
        pdata.push_back(ml->pdata[pdataOffsetNonPadded]);
#else
        int pdataOffsetPadded =
            tools::Vectorizer::SizeOf(ml->nodecount) * i + n;
        int pd = ml->pdata[pdataOffsetPadded];
        int ptype = memb_func[mech->type].dparam_semantics[i];

        // remove extra space added by padding (for pointer to area or ion mech
        // instance)
        if (input_params->branchingDepth > 0 &&
            (ptype == -1 || (ptype > 0 && ptype < 1000))) {
          assert(dataOffsets.at(pd) != -99999);
          pdata.push_back(dataOffsets.at(pd));  // offset to non-padded SoA
                                                // value
        } else
          pdata.push_back(pd);
#endif
      }
      pdata.shrink_to_fit();

      void **vdata = &nt->_vdata[vdataTotalOffset];
      Compartment *compartment = compartments.at(ml->nodeindices[n]);
      assert(compartment->id == ml->nodeindices[n]);

      if (mech->pntMap > 0 || mech->vdataSize > 0)  // vdata
      {
        assert((type == MechanismTypes::kIClamp && mech->vdataSize == 1 &&
                mech->pdataSize == 2 && mech->pntMap > 0) ||
               (type == MechanismTypes::kStochKv && mech->vdataSize == 1 &&
                mech->pdataSize == 5 && mech->pntMap == 0) ||
               ((type == MechanismTypes::kProbAMPANMDA_EMS ||
                 type == MechanismTypes::kProbGABAAB_EMS) &&
                mech->vdataSize == 2 && mech->pdataSize == 3 &&
                mech->pntMap > 0));

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
          const int pointProcOffsetInPdata = 1;
          Point_process *ppn = &nt->pntprocs[pointProcTotalOffset++];
          assert(nt->_vdata[pdata[pointProcOffsetInPdata]] == ppn);
          Point_process *pp = (Point_process *)vdata[0];
          assert(pp != NULL);
          compartment->AddSerializedVdata((unsigned char *)(void *)pp,
                                          sizeof(Point_process));
        }

        if (type == MechanismTypes::kStochKv ||
            type == MechanismTypes::kProbAMPANMDA_EMS ||
            type == MechanismTypes::kProbGABAAB_EMS) {
          int rngOffsetInPdata = mech->type == MechanismTypes::kStochKv ? 3 : 2;
          int rngOffsetInVdata = mech->type == MechanismTypes::kStochKv ? 0 : 1;
          assert(nt->_vdata[pdata[rngOffsetInPdata]] ==
                 vdata[rngOffsetInVdata]);
          nrnran123_State *RNG = (nrnran123_State *)vdata[rngOffsetInVdata];

          // TODO: manual hack: StochKv's current state has NULL pointers, why?
          if (RNG == NULL && type == MechanismTypes::kStochKv)
            compartment->AddSerializedVdata(
                (unsigned char *)new nrnran123_State, sizeof(nrnran123_State));
          else {
            assert(RNG != NULL);
            compartment->AddSerializedVdata((unsigned char *)(void *)RNG,
                                            sizeof(nrnran123_State));
          }
        }
      }
      compartment->AddMechanismInstance(type, n, data.data(), mech->dataSize,
                                        pdata.data(), mech->pdataSize);

      vdataTotalOffset += (unsigned)mech->vdataSize;
    }
    dataTotalOffset += mech->dataSize * ml->nodecount;
    dataTotalPaddedOffset +=
        mech->dataSize * tools::Vectorizer::SizeOf(ml->nodecount);
  }
  for (Compartment *comp : compartments) comp->ShrinkToFit();

  dataOffsets.clear();

  //======= 3 - reconstruct NetCons =====================
  map<neuron_id_t, vector<NetConX *>>
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
    int mechtype = nc->target_->_type;
    int weightscount = pnt_receive_size[mechtype];
    size_t weightsOffset = nc->u.weight_index_;
    assert(weightsOffset < nt->n_weight);

    int nodeid = nrn_threads[nc->target_->_tid]
                     ._ml_list[mechtype]
                     ->nodeindices[nc->target_->_i_instance];
    Compartment *comp = compartments.at(nodeid);
    NetConX *nx = new NetConX(mechtype, (offset_t)nc->target_->_i_instance,
                              (floble_t)nc->delay_, weightsOffset, weightscount,
                              nc->active_);
    double *weights = &nt->weights[weightsOffset];
    comp->AddNetCon(srcgid, nx, weights);
    netcons[srcgid].push_back(nx);
  }

  if (input_params->outputNetconsDot) {
    int netConsFromOthers = 0;
    for (auto nc : netcons) {
      int srcGid = nc.first;
      if (std::find(all_neurons_gids->begin(), all_neurons_gids->end(),
                    srcGid) == all_neurons_gids->end())
        netConsFromOthers++;
      else {
        floble_t minDelay = 99999;
        for (auto ncv : nc.second)  // get minimum delay between neurons
          minDelay = std::min(minDelay, ncv->delay);
        fprintf(fileNetcons, "%d -> %d [label=\"%d (%.2fms)\"];\n", srcGid,
                neuronId, nc.second.size(), minDelay);
      }
    }
#if NEUROX_INPUT_DATALOADER_OUTPUT_EXTERNAL_NETCONS_ == true
    if (netConsFromOthers > 0)
      fprintf(fileNetcons,
              "%s -> %d [label=\"%d\" fontcolor=gray color=gray arrowhead=vee "
              "fontsize=12];\n",
              "external", neuronId, netConsFromOthers);
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
    compartments.at(ppi.nodeId)
        ->AddVecPlay(vec->t_->data(), vec->y_->data(), ppi);
  }

  //======= 5 - recursively create branches tree ===========
  floble_t APthreshold = (floble_t)nrn_threads[nt->id].presyns[0].threshold_;
  int thvar_index = nrn_threads[nt->id].presyns[0].thvar_index_;

  for (Compartment *comp : compartments) comp->ShrinkToFit();

  CreateBranch(nt->id, HPX_NULL, compartments, compartments.at(0),
               ionsInstancesInfo, input_params->branchingDepth, thvar_index,
               APthreshold);

  for (auto c : compartments) delete c;
  for (auto nc : netcons)
    for (auto nvc : nc.second) delete nvc;
  return 0;
}

void DataLoader::CleanCoreneuronData(const bool clean_ion_global_map) {
  nrn_cleanup(clean_ion_global_map);
}

void DataLoader::InitAndLoadCoreneuronData(int argc, char **argv,
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
  NEUROX_MEM_PIN_(uint64_t);

  // To insert mechanisms in the right order, we must first calculate
  // dependencies
  int myNrnThreadsCount = GetMyNrnThreadsCount();

  if (myNrnThreadsCount == 0) NEUROX_MEM_UNPIN_;

  for (int i = 0; i < myNrnThreadsCount; i++) {
    assert(nrn_threads[i].ncell == 1);
  }

  // Different nrn_threads[i] have diff mechanisms; we'll get the union of all
  // neurons' mechs
  std::list<NrnThreadMembList *> orderedMechs;  // list of unique mechanisms
  std::set<int> uniqueMechIds;                  // list of unique mechanism ids

  // insert all mechs from first neuron
  for (NrnThreadMembList *tml = nrn_threads[0].tml; tml != NULL;
       tml = tml->next) {
    orderedMechs.push_back(tml);
    uniqueMechIds.insert(tml->index);
  }

  // insert all mechs from other neurons that do not yet exist, in the right
  // order
  for (int i = 1; i < myNrnThreadsCount; i++)
    for (NrnThreadMembList *tml = nrn_threads[i].tml; tml->next != NULL;
         tml = tml->next)
      if (uniqueMechIds.find(tml->next->index) ==
          uniqueMechIds.end())  // if next mech does not exist
      {  // find correct position in list and insert it there:
        for (auto it = orderedMechs.begin(); it != orderedMechs.end();
             it++)  //...for all existing mechs
          if ((*it)->index == tml->index) {
            auto it_next = std::next(it, 1);
            orderedMechs.insert(it_next, tml->next);  // reminder: .insert adds
                                                      // elements in position
                                                      // before iterator
            // we have it -> it_new -> it_next. Now we will set the value of
            // next pointer
            // auto it_new = std::prev(it_next,1);
            //(*it)->next = *it_new; //Do not change next pointers or neuron is
            // incorrect
            //(*it_new)->next = *it_next;
            uniqueMechIds.insert(tml->next->index);
            break;
          }
      }
  assert(uniqueMechIds.size() == orderedMechs.size());

  std::vector<int> mechIds_serial;
  std::vector<int> dependenciesCount_serial;
  std::vector<int> successorsCount_serial;
  std::vector<int> dependencies_serial;
  std::vector<int> successors_serial;

  for (auto tml_it = orderedMechs.begin(); tml_it != orderedMechs.end();
       tml_it++) {
    auto &tml = *tml_it;
    int type = tml->index;

    vector<int> successors;
    vector<int> dependencies;
    assert(nrn_watch_check[type] == NULL);  // not supported yet

    if (input_params->multiMex) {
      for (auto &tml2 : orderedMechs) {
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
      if (tml->index != orderedMechs.back()->index) {
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
    mechIds_serial.push_back(type);
    dependenciesCount_serial.push_back(dependencies.size());
    successorsCount_serial.push_back(successors.size());
    dependencies_serial.insert(dependencies_serial.end(), dependencies.begin(),
                               dependencies.end());
    successors_serial.insert(successors_serial.end(), successors.begin(),
                             successors.end());
  }

  if (input_params->patternStim[0] != '\0')  //"initialized"
  {
    assert(0);  // not an error: should be initialized already by coreneuron
                // (and above)
    // in the future this should be contidtional (once we get rid of coreneuron
    // data loading)
  }

  // set mechanisms dependencies
  if (neurox::ParallelExecution() && input_params->branchingDepth > 0) {
    // broadcast dependencies, most complete dependency graph will be used
    // across the network
    //(this solves issue of localities loading morphologies without all
    // mechanisms,
    // and processing branches of other localities where those missing
    // mechanisms exist)
    hpx_bcast_rsync(
        DataLoader::SetMechanisms, mechIds_serial.data(),
        sizeof(int) * mechIds_serial.size(), dependenciesCount_serial.data(),
        sizeof(int) * dependenciesCount_serial.size(),
        dependencies_serial.data(), sizeof(int) * dependencies_serial.size(),
        successorsCount_serial.data(),
        sizeof(int) * successorsCount_serial.size(), successors_serial.data(),
        sizeof(int) * successors_serial.size());
  } else {
    // regular setting of mechanisms: all localities have the dependency graph
    // for their morphologies
    DataLoader::SetMechanisms2(
        mechIds_serial.size(), mechIds_serial.data(),
        dependenciesCount_serial.data(), dependencies_serial.data(),
        successorsCount_serial.data(), successors_serial.data());
  }
  NEUROX_MEM_UNPIN_;
}

hpx_action_t DataLoader::Init = 0;
int DataLoader::Init_handler() {
  NEUROX_MEM_PIN_(uint64_t);
  all_neurons_gids = new std::vector<int>();
  all_neurons_mutex = hpx_lco_sema_new(1);

  // even without load balancing, we may require the benchmark info for
  // outputting statistics
  if (hpx_get_my_rank() == 0 &&
      (input_params->loadBalancing || input_params->outputStatistics))
    loadBalancing = new tools::LoadBalancing();

  if (neurox::ParallelExecution()  // disable output of netcons for parallel
                                   // loading
      && input_params->outputNetconsDot) {
    input_params->outputNetconsDot = false;
    if (hpx_get_my_rank() == 0)
      printf("Warning: output of netcons.dot disabled for parallel loading\n");
  }

  if (input_params->outputNetconsDot) {
    assert(HPX_LOCALITIES == 1);
    fileNetcons = fopen(string("netcons.dot").c_str(), "wt");
    fprintf(fileNetcons, "digraph G\n{ bgcolor=%s; layout=circo;\n",
            "transparent");
#if NEUROX_INPUT_DATALOADER_OUTPUT_EXTERNAL_NETCONS_ == true
    fprintf(fileNetcons, "external [color=gray fontcolor=gray];\n");
#endif
  }
  NEUROX_MEM_UNPIN_;
}

hpx_action_t DataLoader::InitNeurons = 0;
int DataLoader::InitNeurons_handler() {
  NEUROX_MEM_PIN_(uint64_t);

  int myNrnThreadsCount = GetMyNrnThreadsCount();

  if (myNrnThreadsCount == 0) NEUROX_MEM_UNPIN_;

#if NEUROX_INPUT_DATALOADER_OUTPUT_CORENEURON_COMPARTMENTS_ == true
  if (inputParams->outputCompartmentsDot) {
    for (int i = 0; i < myNrnThreadsCount; i++) {
      neuron_id_t neuronId = GetNeuronIdFromNrnThreadId(i);
      FILE *fileCompartments =
          fopen(string("compartments_" + to_string(neuronId) + "_NrnThread.dot")
                    .c_str(),
                "wt");
      fprintf(fileCompartments, "graph G%d\n{  node [shape=cylinder];\n",
              neuronId);

      // for all nodes in this NrnThread
      NrnThread *nt = &nrn_threads[i];
      for (int n = nt->ncell; n < nt->end; n++)
        fprintf(fileCompartments, "%d -- %d;\n", nt->_v_parent_index[n], n);
      fprintf(fileCompartments, "}\n");
      fclose(fileCompartments);
    }
  }
#endif

  my_neurons_addr = new std::vector<hpx_t>();
  my_neurons_gids = new std::vector<int>();

  hpx_par_for_sync(DataLoader::CreateNeuron, 0, myNrnThreadsCount, nullptr);

  // all neurons created, advertise them
  hpx_bcast_rsync(DataLoader::AddNeurons, my_neurons_gids->data(),
                  sizeof(int) * my_neurons_gids->size(),
                  my_neurons_addr->data(),
                  sizeof(hpx_t) * my_neurons_addr->size());

  my_neurons_gids->clear();
  delete my_neurons_gids;
  my_neurons_gids = nullptr;
  my_neurons_addr->clear();
  delete my_neurons_addr;
  my_neurons_addr = nullptr;

  if (input_params->allReduceAtLocality) {
    assert(
        0);  // TODO Broken, my_neurons_addrs point to all neurons loaded by me,
    // but can be allocated anywhere
    AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::localityNeurons =
        new std::vector<hpx_t>(my_neurons_addr->begin(),
                               my_neurons_addr->end());
  }

  NEUROX_MEM_UNPIN_;
}

hpx_action_t DataLoader::AddNeurons = 0;
int DataLoader::AddNeurons_handler(const int nargs, const void *args[],
                                   const size_t sizes[]) {
  /**
   * nargs=3 where
   * args[1] = neurons Gids
   * args[2] = neurons hpx addr
   */

  NEUROX_MEM_PIN_(uint64_t);
  assert(nargs == 2);
  assert(sizes[0] / sizeof(int) == sizes[1] / sizeof(hpx_t));

  const int recv_neurons_count = sizes[0] / sizeof(int);
  const int *neurons_gids = (const int *)args[0];
  const hpx_t *neurons_addr = (const hpx_t *)args[1];

  hpx_lco_sema_p(all_neurons_mutex);

  all_neurons_gids->insert(all_neurons_gids->end(), neurons_gids,
                           neurons_gids + recv_neurons_count);

  // add the values to neurox::neurons
  hpx_t *neurons_new = new hpx_t[neurox::neurons_count + recv_neurons_count];
  copy(neurox::neurons, neurox::neurons + neurox::neurons_count, neurons_new);
  copy(neurons_addr, neurons_addr + recv_neurons_count,
       neurons_new + neurox::neurons_count);

  delete[] neurox::neurons;
  neurox::neurons_count += recv_neurons_count;
  neurox::neurons = neurons_new;

  assert(all_neurons_gids->size() == neurox::neurons_count);
  hpx_lco_sema_v_sync(all_neurons_mutex);
  NEUROX_MEM_UNPIN_;
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

  NEUROX_MEM_PIN_(uint64_t);
  assert(nargs == 5);

  const int *orderedMechs = (const int *)args[0];
  const int *dependenciesCount = (const int *)args[1];
  const int *dependencies = (const int *)args[2];
  const int *successorsCount = (const int *)args[3];
  const int *successors = (const int *)args[4];

  const int mechsCount = sizes[0] / sizeof(int);

  SetMechanisms2(mechsCount, orderedMechs, dependenciesCount, dependencies,
                 successorsCount, successors);

  NEUROX_MEM_UNPIN_;
}

void DataLoader::SetMechanisms2(const int mechsCount, const int *mechIds,
                                const int *dependenciesCount,
                                const int *dependencies,
                                const int *successorsCount,
                                const int *successors) {
  hpx_lco_sema_p(all_neurons_mutex);

  // if I'm receiving new information, rebuild mechanisms
  if (mechsCount > neurox::mechanisms_count) {
    // delete existing data (if any)
    if (neurox::mechanisms_count > 0) {
      for (int m = 0; m < neurox::mechanisms_count; m++) delete mechanisms[m];
      delete[] neurox::mechanisms;
      delete[] neurox::mechanisms_map;
    }

    neurox::mechanisms_count = mechsCount;
    neurox::mechanisms_map = new int[n_memb_func];
    neurox::mechanisms = new Mechanism *[mechsCount];

    for (int i = 0; i < n_memb_func; i++) neurox::mechanisms_map[i] = -1;

    int dependenciesOffset = 0;
    int successorsOffset = 0;
    for (int m = 0; m < mechsCount; m++) {
      int type = mechIds[m];
      neurox::mechanisms_map[type] = m;

      int symLength =
          memb_func[type].sym ? std::strlen(memb_func[type].sym) : 0;
      const int *mechDependencies = dependenciesCount[m] == 0
                                        ? nullptr
                                        : &dependencies[dependenciesOffset];
      const int *mechSuccessors =
          successorsCount[m] == 0 ? nullptr : &successors[successorsOffset];
      Mechanism *mech = new Mechanism(
          type, nrn_prop_param_size_[type], nrn_prop_dparam_size_[type],
          nrn_is_artificial_[type], pnt_map[type], nrn_is_ion(type), symLength,
          memb_func[type].sym, memb_func[type], dependenciesCount[m],
          mechDependencies, successorsCount[m], mechSuccessors);
      neurox::mechanisms[m] = mech;

      successorsOffset += successorsCount[m];
      dependenciesOffset += dependenciesCount[m];
    }

    // set parent ion index
    for (int m = 0; m < mechsCount; m++) {
      if (input_params->multiMex) {
        Mechanism *mech = neurox::mechanisms[m];
        for (int d = 0; d < mech->dependenciesCount; d++) {
          Mechanism *parent = GetMechanismFromType(mech->dependencies[d]);
          if (strcmp("SK_E2", mech->membFunc.sym) == 0 &&
              strcmp("ca_ion", parent->membFunc.sym) == 0)
            continue;  // TODO hard coded exception

          if (parent->GetIonIndex() < Mechanism::IonTypes::kSizeWriteableIons)
            mech->dependencyIonIndex = parent->GetIonIndex();
        }
      }
    }
  }
  hpx_lco_sema_v_sync(all_neurons_mutex);
}

hpx_action_t DataLoader::Finalize = 0;
int DataLoader::Finalize_handler() {
  NEUROX_MEM_PIN_(uint64_t);

  if (input_params->allReduceAtLocality)
    AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::localityNeurons
        ->clear();
  delete AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::localityNeurons;

  if (input_params->outputNetconsDot) {
    fprintf(fileNetcons, "}\n");
    fclose(fileNetcons);
  }

  if (input_params->outputMechanismsDot) {
    FILE *fileMechs =
        fopen(string("mechanisms_" + std::to_string(hpx_get_my_rank()) + ".dot")
                  .c_str(),
              "wt");
    fprintf(fileMechs, "digraph G\n{ bgcolor=%s; %s\n", "transparent",
            !input_params->multiMex ? "layout=circo; scale=0.23;" : "");
    fprintf(fileMechs, "graph [ratio=0.3];\n", "start");
    fprintf(fileMechs, "%s [style=filled, shape=Mdiamond, fillcolor=beige];\n",
            "start");
    fprintf(fileMechs, "%s [style=filled, shape=Mdiamond, fillcolor=beige];\n",
            "end");
    fprintf(fileMechs, "\"%s (%d)\" [style=filled, fillcolor=beige];\n",
            GetMechanismFromType(CAP)->membFunc.sym, CAP);
    if (!input_params->multiMex) {
      fprintf(fileMechs, "end -> start [color=transparent];\n");
      fprintf(fileMechs, "start -> \"%s (%d)\";\n",
              GetMechanismFromType(CAP)->membFunc.sym, CAP);
    }
    for (int m = 0; m < mechanisms_count; m++) {
      Mechanism *mech = neurox::mechanisms[m];

      if (mech->pntMap > 0)  // if is point process make it dotted
        fprintf(fileMechs, "\"%s (%d)\" [style=dashed];\n", mech->membFunc.sym,
                mech->type);

      if (mech->dependenciesCount == 0 && mech->type != CAP)  // top mechanism
        fprintf(fileMechs, "%s -> \"%s (%d)\";\n", "start", mech->membFunc.sym,
                mech->type);

      if (mech->successorsCount == 0 && mech->type != CAP)  // bottom mechanism
        fprintf(fileMechs, "\"%s (%d)\" -> %s;\n", mech->membFunc.sym,
                mech->type, "end");

      for (int s = 0; s < mech->successorsCount; s++) {
        Mechanism *successor = GetMechanismFromType(mech->successors[s]);
        fprintf(fileMechs, "\"%s (%d)\" -> \"%s (%d)\";\n", mech->membFunc.sym,
                mech->type, successor->membFunc.sym, successor->type);
      }

      if (input_params->multiMex)
        for (int d = 0; d < mech->dependenciesCount; d++) {
          Mechanism *parent = GetMechanismFromType(mech->dependencies[d]);
          if (strcmp("SK_E2", mech->membFunc.sym) == 0 &&
              strcmp("ca_ion", parent->membFunc.sym) == 0)
            continue;  // TODO: hardcoded exception
          if (parent->GetIonIndex() <
              Mechanism::IonTypes::kSizeWriteableIons)  // ie is writeable
            fprintf(
                fileMechs,
                "\"%s (%d)\" -> \"%s (%d)\" [style=dashed, arrowtype=open];\n",
                mech->membFunc.sym, mech->type, parent->membFunc.sym,
                parent->type);
        }
    }
    fprintf(fileMechs, "}\n");
    fclose(fileMechs);
  }

#ifndef NDEBUG
  if (HPX_LOCALITY_ID == 0) {
    for (int m = 0; m < neurox::mechanisms_count; m++) {
      Mechanism *mech = neurox::mechanisms[m];
      printf(
          "- %s (%d), dataSize %d, pdataSize %d, isArtificial %d, pntMap %d, "
          "isIon %d, symLength %d, %d successors, %d dependencies\n",
          mech->membFunc.sym, mech->type, mech->dataSize, mech->pdataSize,
          mech->isArtificial, mech->pntMap, mech->isIon, mech->symLength,
          mech->successorsCount, mech->dependenciesCount);
    }
  }
#endif

  all_neurons_gids->clear();
  delete all_neurons_gids;
  all_neurons_gids = nullptr;

  nrn_setup_cleanup();
  hpx_lco_delete_sync(all_neurons_mutex);

#if defined(NDEBUG)
  // if not on debug, there's no CoreNeuron comparison, so data can be
  // cleaned-up now, except global_ion_map
  DataLoader::CleanCoreneuronData(false);
#endif

  if (input_params->outputStatistics)
    tools::LoadBalancing::LoadBalancing::PrintTable();

  delete loadBalancing;
  loadBalancing = nullptr;
  NEUROX_MEM_UNPIN_;
}

void DataLoader::GetAllChildrenCompartments(deque<Compartment *> &subSection,
                                            Compartment *topCompartment) {
  // parent added first, to respect solvers logic, of parents' ids first
  subSection.push_back(topCompartment);
  for (int c = 0; c < topCompartment->branches.size(); c++)
    GetAllChildrenCompartments(subSection, topCompartment->branches.at(c));
}

void DataLoader::GetMechInstanceMap(deque<Compartment *> &compartments,
                                    vector<map<int, int>> &mechsInstancesMap) {
  vector<deque<int>> mechsInstancesIds(
      neurox::mechanisms_count);  // mech offset -> list of mech instance id
  for (Compartment *comp : compartments)
    for (int m = 0; m < comp->mechsTypes.size(); m++)  // for all instances
    {
      int type = comp->mechsTypes.at(m);
      int mechOffset = neurox::mechanisms_map[type];
      mechsInstancesIds[mechOffset].push_back(comp->mechsInstances[m]);
    }

  // convert neuron mech-instances ids from neuron- to branch-level
  for (int m = 0; m < neurox::mechanisms_count; m++)
    for (int i = 0; i < mechsInstancesIds[m].size(); i++) {
      int oldInstanceId = mechsInstancesIds.at(m).at(i);
      mechsInstancesMap[m][oldInstanceId] = i;
    }
}

void DataLoader::GetNetConsBranchData(deque<Compartment *> &compartments,
                                      vector<NetConX> &branchNetCons,
                                      vector<neuron_id_t> &branchNetConsPreId,
                                      vector<floble_t> &branchWeights,
                                      vector<map<int, int>> *mechInstanceMap) {
  for (auto &comp : compartments) {
    branchNetCons.insert(branchNetCons.end(), comp->netcons.begin(),
                         comp->netcons.end());
    branchNetConsPreId.insert(branchNetConsPreId.end(),
                              comp->netconsPreSynIds.begin(),
                              comp->netconsPreSynIds.end());
    branchWeights.insert(branchWeights.end(), comp->netconsWeights.begin(),
                         comp->netconsWeights.end());
  }

  // correct weighIndex variables
  int weightOffset = 0;
  for (NetConX &nc : branchNetCons) {
    nc.weightIndex = weightOffset;
    weightOffset += nc.weightsCount;
  }

  // convert mech instance id from neuron to branch level
  if (mechInstanceMap)
    for (NetConX &nc : branchNetCons)
      nc.mechInstance = (*mechInstanceMap)[neurox::mechanisms_map[nc.mechType]]
                                          [nc.mechInstance];
}

void DataLoader::GetVecPlayBranchData(deque<Compartment *> &compartments,
                                      vector<floble_t> &vecPlayTdata,
                                      vector<floble_t> &vecPlayYdata,
                                      vector<PointProcInfo> &vecPlayInfo,
                                      vector<map<int, int>> *mechInstanceMap) {
  // convert node id and mech instance id in PointProcess from neuron to branch
  // level
  if (mechInstanceMap) {
    std::map<int, int> fromOldToNewCompartmentId;
    for (int n = 0; n < compartments.size(); n++)
      fromOldToNewCompartmentId[compartments.at(n)->id] = n;

    for (int p = 0; p < vecPlayInfo.size(); p++) {
      PointProcInfo &ppi = vecPlayInfo[p];
      ppi.mechInstance =
          (offset_t)(*mechInstanceMap)[neurox::mechanisms_map[ppi.mechType]]
                                      [ppi.mechInstance];
      ppi.nodeId = fromOldToNewCompartmentId[ppi.nodeId];
    }
  }

  for (auto comp : compartments) {
    vecPlayTdata.insert(vecPlayTdata.end(), comp->vecPlayTdata.begin(),
                        comp->vecPlayTdata.end());
    vecPlayYdata.insert(vecPlayYdata.end(), comp->vecPlayYdata.begin(),
                        comp->vecPlayYdata.end());
    vecPlayInfo.insert(vecPlayInfo.end(), comp->vecPlayInfo.begin(),
                       comp->vecPlayInfo.end());
  }
}

int DataLoader::GetBranchData(
    deque<Compartment *> &compartments, vector<floble_t> &data,
    vector<offset_t> &pdata, vector<unsigned char> &vdata, vector<offset_t> &p,
    vector<offset_t> &instancesCount, vector<offset_t> &nodesIndices, int N,
    vector<DataLoader::IonInstancesInfo> &ionsInstancesInfo,
    vector<map<int, int>> *mechInstanceMap) {
  for (auto comp : compartments) {
    assert(comp != NULL);
  }

  int n = 0;  // number of compartments
  int vdataPointerOffset = 0;

  ////// Basic information for RHS, D, A, B, V and area
  for (auto comp : compartments) data.push_back(comp->rhs);
  for (auto comp : compartments) data.push_back(comp->d);
  for (auto comp : compartments) data.push_back(comp->a);
  for (auto comp : compartments) data.push_back(comp->b);
  for (auto comp : compartments) data.push_back(comp->v);
  for (auto comp : compartments) data.push_back(comp->area);
  for (auto comp : compartments) p.push_back(comp->p);

  ////// Tree of neurons: convert from neuron- to branch-level
  std::map<int, int> fromOldToNewCompartmentId;
  for (Compartment *comp : compartments) {
    fromOldToNewCompartmentId[comp->id] = n;
    comp->id = n;
    n++;
  }
  if (mechInstanceMap) {
    p.at(0) = 0;  // top node gets parent Id 0 as in Coreneuron
    for (int i = 1; i < p.size(); i++)
      p.at(i) = fromOldToNewCompartmentId.at(p.at(i));
  }

  ////// Mechanisms instances: merge all instances of all compartments into
  /// instances of the branch
  vector<vector<floble_t>> dataMechs(neurox::mechanisms_count);
  vector<vector<offset_t>> pdataMechs(neurox::mechanisms_count);
  vector<vector<offset_t>> nodesIndicesMechs(neurox::mechanisms_count);
  vector<vector<unsigned char>> vdataMechs(neurox::mechanisms_count);
  vector<deque<int>> mechsInstancesIds(
      neurox::mechanisms_count);  // mech-offset -> list of pointers to mech
                                  // instance value

  map<pair<int, offset_t>, offset_t> ionInstanceToDataOffset;  // from pair of <
                                                               // ion mech type,
                                                               // OLD node id>
                                                               // to ion offset
                                                               // in NEW
                                                               // representation

  for (Compartment *comp : compartments) {
    int compDataOffset = 0;
    int compPdataOffset = 0;
    int compVdataOffset = 0;
    for (int m = 0; m < comp->mechsTypes.size(); m++)  // for all instances
    {
      int type = comp->mechsTypes[m];
      int mechOffset = mechanisms_map[type];
      assert(mechOffset >= 0 && mechOffset < mechanisms_count);
      Mechanism *mech = mechanisms[mechOffset];
      dataMechs[mechOffset].insert(
          dataMechs[mechOffset].end(), &comp->data[compDataOffset],
          &comp->data[compDataOffset + mech->dataSize]);
      pdataMechs[mechOffset].insert(
          pdataMechs[mechOffset].end(), &comp->pdata[compPdataOffset],
          &comp->pdata[compPdataOffset + mech->pdataSize]);
      nodesIndicesMechs[mechOffset].push_back(comp->id);
      instancesCount[mechOffset]++;
      compDataOffset += mech->dataSize;
      compPdataOffset += mech->pdataSize;
      if (mech->pntMap > 0 || mech->vdataSize > 0) {
        assert((type == MechanismTypes::kIClamp && mech->vdataSize == 1 &&
                mech->pdataSize == 2 && mech->pntMap > 0) ||
               (type == MechanismTypes::kStochKv && mech->vdataSize == 1 &&
                mech->pdataSize == 5 && mech->pntMap == 0) ||
               ((type == MechanismTypes::kProbAMPANMDA_EMS ||
                 type == MechanismTypes::kProbGABAAB_EMS) &&
                mech->vdataSize == 2 && mech->pdataSize == 3 &&
                mech->pntMap > 0));

        size_t totalVdataSize = 0;
        if (type == MechanismTypes::kIClamp ||
            type == MechanismTypes::kProbAMPANMDA_EMS ||
            type == MechanismTypes::kProbGABAAB_EMS)
          totalVdataSize += sizeof(Point_process);
        if (type == MechanismTypes::kStochKv ||
            type == MechanismTypes::kProbAMPANMDA_EMS ||
            type == MechanismTypes::kProbGABAAB_EMS)
          totalVdataSize += sizeof(nrnran123_State);
        vdataMechs[mechOffset].insert(
            vdataMechs[mechOffset].end(), &comp->vdata[compVdataOffset],
            &comp->vdata[compVdataOffset + totalVdataSize]);
        compVdataOffset += totalVdataSize;
      }
    }
  }

  // merge all mechanisms vectors in the final one
  // store the offset of each mechanism data (for later)
  for (int m = 0; m < neurox::mechanisms_count; m++) {
    Mechanism *mech = neurox::mechanisms[m];
    int dataOffset = 0;
    int pdataOffset = 0;
    int vdataOffset = 0;
    for (int i = 0; i < nodesIndicesMechs[m].size(); i++)  // for all instances
    {
      assert(ionInstanceToDataOffset.find(
                 make_pair(mech->type, nodesIndicesMechs[m][i])) ==
             ionInstanceToDataOffset.end());
      if (mech->isIon &&
          input_params->branchingDepth > 0)  // for pdata calculation
        ionInstanceToDataOffset[make_pair(
            mech->type, nodesIndicesMechs[m][i])] = data.size();
      data.insert(data.end(), &dataMechs[m][dataOffset],
                  &dataMechs[m][dataOffset + mech->dataSize]);
      nodesIndices.push_back(nodesIndicesMechs[m][i]);
      dataOffset += mech->dataSize;

      if (mech->pntMap > 0 || mech->vdataSize > 0) {
        int totalVdataSize = 0;
        if (mech->type == MechanismTypes::kIClamp ||
            mech->type == MechanismTypes::kProbAMPANMDA_EMS ||
            mech->type == MechanismTypes::kProbGABAAB_EMS)
          totalVdataSize += sizeof(Point_process);
        if (mech->type == MechanismTypes::kStochKv ||
            mech->type == MechanismTypes::kProbAMPANMDA_EMS ||
            mech->type == MechanismTypes::kProbGABAAB_EMS)
          totalVdataSize += sizeof(nrnran123_State);
        assert(totalVdataSize > 0);
        vdata.insert(vdata.end(), &vdataMechs[m][vdataOffset],
                     &vdataMechs[m][vdataOffset + totalVdataSize]);
        vdataOffset += totalVdataSize;
      }

      if (input_params->branchingDepth >
          0)  // if we need to recalculate offsets or remove padding
      {
        for (int p = pdataOffset; p < pdataOffset + mech->pdataSize; p++) {
          offset_t pd = pdataMechs.at(m).at(p);
          int ptype = memb_func[mech->type].dparam_semantics[p - pdataOffset];
          switch (ptype) {
            case -1:  //"area" (6th field)
            {
              assert(pd >= N * 5 && pd < N * 6);
              offset_t oldId = pd - N * 5;
              offset_t newId = fromOldToNewCompartmentId.at(oldId);
              pdataMechs.at(m).at(p) = n * 5 + newId;
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
              pdataMechs.at(m).at(p) = (offset_t)vdataPointerOffset++;
              break;
            case -8:      //"bbcorepointer"
              assert(0);  // watch condition, not supported
              break;
            default:
              if (ptype > 0 && ptype < 1000)  // name preffixed by '#'
              {
                // ptype is the ion (mechanism) type it depends on
                // pdata is an offset in nt->data (a var in the ion)
                // pd points to SoA notation, independently of the LAYOUT
                // (converted before)

                Mechanism *ion = neurox::GetMechanismFromType(ptype);
                Mechanism::IonTypes ionOffset = ion->GetIonIndex();
                IonInstancesInfo &ionInfo =
                    ionsInstancesInfo.at((int)ionOffset);
                int dataStart = ionInfo.dataStart;
                assert(pd >= dataStart && pd < ionInfo.dataEnd);

                int instanceOffset =
                    trunc((double)(pd - dataStart) / (double)ion->dataSize);
                int instanceVariableOffset = (pd - dataStart) % ion->dataSize;
                int nodeId = ionInfo.nodeIds.at(instanceOffset);
                int newNodeId = fromOldToNewCompartmentId.at(nodeId);
                pdataMechs.at(m).at(p) = ionInstanceToDataOffset.at(
                                             make_pair(ion->type, newNodeId)) +
                                         instanceVariableOffset;
                assert(pdataMechs.at(m).at(p) >= n * 6);
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
      pdata.insert(pdata.end(), &pdataMechs[m][pdataOffset],
                   &pdataMechs[m][pdataOffset + mech->pdataSize]);
      pdataOffset += mech->pdataSize;
    }
    dataMechs[m].clear();
    pdataMechs[m].clear();
    vdataMechs[m].clear();
    nodesIndicesMechs[m].clear();

    // convert neuron mech-instances ids from neuron- to branch-level
    if (mechInstanceMap)
      for (int i = 0; i < mechsInstancesIds[m].size(); i++) {
        int oldInstanceId = mechsInstancesIds[m][i];
        (*mechInstanceMap)[m][oldInstanceId] = i;
      }
  }
  return n;
}

bool CompareCompartmentsPtrsIds(Compartment *a, Compartment *b) {
  return a->id < b->id;
}

hpx_t DataLoader::CreateBranch(
    int nrnThreadId, hpx_t somaBranchAddr,
    deque<Compartment *> &allCompartments, Compartment *topCompartment,
    vector<DataLoader::IonInstancesInfo> &ionsInstancesInfo, int branchingDepth,
    int thvar_index /*AIS*/, floble_t APthreshold /*AIS*/) {
  assert(topCompartment != NULL);
  offset_t n;              // number of compartments in branch
  vector<floble_t> data;   // compartments info (RHS, D, A, B, V, AREA)*n
  vector<offset_t> pdata;  // pointers to data
  vector<offset_t> p;      // parent nodes index
  vector<offset_t> instancesCount(mechanisms_count);
  vector<offset_t> nodesIndices;
  vector<hpx_t> branches;
  vector<unsigned char> vdata;  // Serialized Point Processes and Random123
  int N = allCompartments.size();

  bool isSoma = somaBranchAddr == HPX_NULL;

  // Vector Play instances
  vector<floble_t> vecPlayT;
  vector<floble_t> vecPlayY;
  vector<PointProcInfo> vecPlayInfo;

  // branch NetCons
  vector<NetConX> branchNetCons;
  vector<neuron_id_t> branchNetConsPreId;
  vector<floble_t> branchWeights;

  Compartment *bottomCompartment = nullptr;  // iterator for subsections

  if (input_params->branchingDepth == 0)  // Flat a la Coreneuron
  {
    n = GetBranchData(allCompartments, data, pdata, vdata, p, instancesCount,
                      nodesIndices, N, ionsInstancesInfo, NULL);
    GetVecPlayBranchData(allCompartments, vecPlayT, vecPlayY, vecPlayInfo,
                         NULL);
    GetNetConsBranchData(allCompartments, branchNetCons, branchNetConsPreId,
                         branchWeights, NULL);
  } else if (input_params->branchingDepth > 0)  // branch-parallelism
  {
    deque<Compartment *> subSection;
    if (branchingDepth > 0)  // node of a tree, with branches
    {
      // subsection is the set of all sequential compartments until finding a
      // bifurcation
      subSection.push_back(topCompartment);
      for (bottomCompartment = topCompartment;
           bottomCompartment->branches.size() == 1;
           bottomCompartment = bottomCompartment->branches.front())
        subSection.push_back(bottomCompartment->branches.front());
    } else  // leaf of the tree
    {
      // subsection is the set of all children compartments (recursively)
      GetAllChildrenCompartments(subSection, topCompartment);
    }

    // this step is only necesary so that data has the same alignment as
    // CoreNeuron
    // and allows one to compare results . (can be removed for non-debug mode)
    std::sort(subSection.begin(), subSection.end(), CompareCompartmentsPtrsIds);

    // create sub-section of branch
    vector<map<int, int>> mechInstanceMap(mechanisms_count);  // mech-offset ->
                                                              // ( map[old
                                                              // instance]->to
                                                              // new instance )
    GetMechInstanceMap(subSection, mechInstanceMap);
    n = GetBranchData(subSection, data, pdata, vdata, p, instancesCount,
                      nodesIndices, N, ionsInstancesInfo, &mechInstanceMap);
    GetVecPlayBranchData(subSection, vecPlayT, vecPlayY, vecPlayInfo,
                         &mechInstanceMap);
    GetNetConsBranchData(subSection, branchNetCons, branchNetConsPreId,
                         branchWeights, &mechInstanceMap);
  }

  int neuronRank = hpx_get_my_rank();
  if (input_params->loadBalancing || input_params->outputStatistics) {
    // Benchmark and assign this branch to least busy compute node (except soma
    // and AIS)
    // Note: we do this after children creation so that we use top (lighter)
    // branches to balance work load
    hpx_t tempBranchAddr =
        hpx_gas_alloc_local(1, sizeof(Branch), NEUROX_MEM_ALIGNMENT_);
    bool runBenchmarkAndClear = true;
    int dumbThresholdOffset = 0;
    double timeElapsed = -1;
    hpx_call_sync(
        tempBranchAddr, Branch::Init, &timeElapsed,
        sizeof(timeElapsed),  // output
        &n, sizeof(offset_t), &nrnThreadId, sizeof(int), &dumbThresholdOffset,
        sizeof(int), data.size() > 0 ? data.data() : nullptr,
        sizeof(floble_t) * data.size(),
        pdata.size() > 0 ? pdata.data() : nullptr,
        sizeof(offset_t) * pdata.size(), instancesCount.data(),
        instancesCount.size() * sizeof(offset_t), nodesIndices.data(),
        nodesIndices.size() * sizeof(offset_t), &somaBranchAddr, sizeof(hpx_t),
        nullptr, 0,                             // no branches
        p.data(), sizeof(offset_t) * p.size(),  // force use of parent index
        vecPlayT.size() > 0 ? vecPlayT.data() : nullptr,
        sizeof(floble_t) * vecPlayT.size(),
        vecPlayY.size() > 0 ? vecPlayY.data() : nullptr,
        sizeof(floble_t) * vecPlayY.size(),
        vecPlayInfo.size() > 0 ? vecPlayInfo.data() : nullptr,
        sizeof(PointProcInfo) * vecPlayInfo.size(),
        branchNetCons.size() > 0 ? branchNetCons.data() : nullptr,
        sizeof(NetConX) * branchNetCons.size(),
        branchNetConsPreId.size() > 0 ? branchNetConsPreId.data() : nullptr,
        sizeof(neuron_id_t) * branchNetConsPreId.size(),
        branchWeights.size() > 0 ? branchWeights.data() : nullptr,
        sizeof(floble_t) * branchWeights.size(),
        vdata.size() > 0 ? vdata.data() : nullptr,
        sizeof(unsigned char) * vdata.size(), &runBenchmarkAndClear,
        sizeof(bool));
    assert(timeElapsed > 0);
    hpx_gas_clear_affinity(tempBranchAddr);

    if (input_params->loadBalancing) {
      // ask master rank to query load balancing table and tell me where to
      // allocate this branch
      hpx_call_sync(HPX_THERE(0), tools::LoadBalancing::QueryLoadBalancingTable,
                    &neuronRank, sizeof(int),       // output
                    &timeElapsed, sizeof(double));  // input [0]
    } else if (input_params->outputStatistics)      // no load balancing
    {
      // compute node already decided, tell master rank to store benchmark info,
      // to be able to output statistics
      hpx_call_sync(HPX_THERE(0), tools::LoadBalancing::QueryLoadBalancingTable,
                    nullptr, 0,                    // output
                    &timeElapsed, sizeof(double),  // input[0]
                    &neuronRank, sizeof(int));     // input[1]
    }

#ifndef NDEBUG
    printf("- %s of neuron nrn_id %d allocated to rank %d (%.6f ms)\n",
           isSoma ? "soma" : (thvar_index != -1 ? "AIS" : "dendrite"),
           nrnThreadId, neuronRank, timeElapsed);
#endif
  }

  // allocate and initialize branch on the respective owner
  hpx_t branchAddr = hpx_gas_alloc_local_at_sync(
      1, sizeof(Branch), NEUROX_MEM_ALIGNMENT_, HPX_THERE(neuronRank));

  // update hpx address of soma
  somaBranchAddr = isSoma ? branchAddr : somaBranchAddr;

  // allocate children branches recursively (if any)
  if (bottomCompartment)
    for (size_t c = 0; c < bottomCompartment->branches.size(); c++)
      branches.push_back(CreateBranch(
          nrnThreadId, somaBranchAddr, allCompartments,
          bottomCompartment->branches[c], ionsInstancesInfo, branchingDepth - 1,
          isSoma && c == 0 ? thvar_index - n
                           : -1)); /*offset in AIS = offset in soma - nt->end */

  // if branching, soma has not threshold var (it was past to AIS above)
  thvar_index = isSoma && input_params->branchingDepth > 0 ? -1 : thvar_index;

  hpx_call_sync(
      branchAddr, Branch::Init, NULL, 0,  // no timing
      &n, sizeof(offset_t), &nrnThreadId, sizeof(int), &thvar_index,
      sizeof(int), data.size() > 0 ? data.data() : nullptr,
      sizeof(floble_t) * data.size(), pdata.size() > 0 ? pdata.data() : nullptr,
      sizeof(offset_t) * pdata.size(), instancesCount.data(),
      instancesCount.size() * sizeof(offset_t), nodesIndices.data(),
      nodesIndices.size() * sizeof(offset_t), &somaBranchAddr, sizeof(hpx_t),
      branches.size() ? branches.data() : nullptr,
      branches.size() ? sizeof(hpx_t) * branches.size() : 0,
      branches.size() ? nullptr : p.data(),
      branches.size() ? 0 : sizeof(offset_t) * p.size(),
      vecPlayT.size() > 0 ? vecPlayT.data() : nullptr,
      sizeof(floble_t) * vecPlayT.size(),
      vecPlayY.size() > 0 ? vecPlayY.data() : nullptr,
      sizeof(floble_t) * vecPlayY.size(),
      vecPlayInfo.size() > 0 ? vecPlayInfo.data() : nullptr,
      sizeof(PointProcInfo) * vecPlayInfo.size(),
      branchNetCons.size() > 0 ? branchNetCons.data() : nullptr,
      sizeof(NetConX) * branchNetCons.size(),
      branchNetConsPreId.size() > 0 ? branchNetConsPreId.data() : nullptr,
      sizeof(neuron_id_t) * branchNetConsPreId.size(),
      branchWeights.size() > 0 ? branchWeights.data() : nullptr,
      sizeof(floble_t) * branchWeights.size(),
      vdata.size() > 0 ? vdata.data() : nullptr,
      sizeof(unsigned char) * vdata.size());

  if (isSoma) {
    // create soma data structure
    int neuronId = GetNeuronIdFromNrnThreadId(nrnThreadId);
    hpx_call_sync(branchAddr, Branch::InitSoma, NULL, 0, &neuronId,
                  sizeof(neuron_id_t), &APthreshold, sizeof(floble_t));

    hpx_lco_sema_p(all_neurons_mutex);
    my_neurons_addr->push_back(branchAddr);
    my_neurons_gids->push_back(neuronId);
    hpx_lco_sema_v_sync(all_neurons_mutex);
  }
  return branchAddr;
}

hpx_action_t DataLoader::InitNetcons = 0;
int DataLoader::InitNetcons_handler() {
  NEUROX_MEM_PIN_(Branch);
  NEUROX_RECURSIVE_BRANCH_ASYNC_CALL_(DataLoader::InitNetcons);

  if (local->soma && input_params->outputNetconsDot)
    fprintf(fileNetcons, "%d [style=filled, shape=ellipse];\n",
            local->soma->gid);

  std::deque<std::pair<hpx_t, floble_t>>
      netcons;  // set of <srcAddr, minDelay> synapses to notify
  std::deque<std::pair<int, spike_time_t>>
      dependencies;  // set of <srcGid, nextNotificationtime> for dependencies

  const floble_t impossiblyLargeDelay = 99999999;
  for (int i = 0; i < all_neurons_gids->size();
       i++)  // loop through all neurons
  {
    neuron_id_t srcGid = all_neurons_gids->at(i);

    // if I'm connected to it (ie is not artificial or non-existent)
    if (local->netcons.find(srcGid) != local->netcons.end()) {
      hpx_t srcAddr = neurox::neurons[i];
      floble_t minDelay = impossiblyLargeDelay;
      for (NetConX *&nc : local->netcons.at(srcGid))
        if (nc->active) minDelay = min(minDelay, nc->delay);

#if NETCONS_OUTPUT_ADDITIONAL_VALIDATION_FILE == true
      if (input_params->outputNetconsDot && minDelay != impossiblyLargeDelay)
        fprintf(fileNetcons, "%d -> %d [label=\"%d (%.2fms)\"];\n", srcGid,
                local->soma->gid, local->netcons.at(srcGid).size(), minDelay);
#endif

      // tell the neuron to add a synapse to this branch and inform him of the
      // fastest netcon we have
      if (minDelay != impossiblyLargeDelay)  // if any active synapse
      {
        // add this netcon to list of netcons to be communicated
        netcons.push_back(make_pair(srcAddr, minDelay));

        // add this pre-syn neuron as my time-dependency
        if (input_params->algorithm == AlgorithmType::kBenchmarkAll ||
            input_params->algorithm ==
                AlgorithmType::kBackwardEulerTimeDependencyLCO) {
          spike_time_t notificationTime =
              input_params->tstart +
              minDelay * TimeDependencyLCOAlgorithm::TimeDependencies::
                             notificationIntervalRatio;
          dependencies.push_back(make_pair(srcGid, notificationTime));
        }
      }
    }
  }

  // inform pre-syn neuron that he connects to me
  hpx_t netconsLCO = hpx_lco_and_new(netcons.size());
  hpx_t topBranchAddr = local->soma ? target : local->branchTree->topBranchAddr;
  int myGid = local->soma ? local->soma->gid : -1;
  for (std::pair<hpx_t, floble_t> &nc : netcons)
    hpx_call(nc.first, DataLoader::AddSynapse, netconsLCO, &target,
             sizeof(hpx_t), &nc.second, sizeof(nc.second), &topBranchAddr,
             sizeof(hpx_t), &myGid, sizeof(int));
  hpx_lco_wait_reset(netconsLCO);
  hpx_lco_delete_sync(netconsLCO);

  // inform my soma of my time dependencies
  hpx_t dependenciesLCO = hpx_lco_and_new(dependencies.size());
  bool initPhase = true;
  for (std::pair<int, spike_time_t> &dep : dependencies)
    hpx_call(topBranchAddr, Branch::UpdateTimeDependency, dependenciesLCO,
             &dep.first, sizeof(neuron_id_t), &dep.second, sizeof(spike_time_t),
             &initPhase, sizeof(bool));
  hpx_lco_wait_reset(dependenciesLCO);
  hpx_lco_delete_sync(dependenciesLCO);

  NEUROX_RECURSIVE_BRANCH_ASYNC_WAIT_;
  NEUROX_MEM_UNPIN_;
}

hpx_action_t DataLoader::AddSynapse = 0;
int DataLoader::AddSynapse_handler(const int nargs, const void *args[],
                                   const size_t[]) {
  NEUROX_MEM_PIN_(Branch);
  assert(local->soma);
  assert(nargs == 4);
  hpx_t addr = *(const hpx_t *)args[0];
  floble_t minDelay = *(const floble_t *)args[1];
  hpx_t topBranchAddr = *(const hpx_t *)args[2];
  int destinationGid = *(const int *)args[3];
  local->soma->AddSynapse(
      new Neuron::Synapse(addr, minDelay, topBranchAddr, destinationGid));
  NEUROX_MEM_UNPIN_;
}

void DataLoader::RegisterHpxActions() {
  NEUROX_REGISTER_ACTION_(NEUROX_ACTION_ZERO_VAR_, DataLoader::Init);
  NEUROX_REGISTER_ACTION_(NEUROX_ACTION_ZERO_VAR_, DataLoader::InitMechanisms);
  NEUROX_REGISTER_ACTION_(NEUROX_ACTION_ZERO_VAR_, DataLoader::InitNeurons);
  NEUROX_REGISTER_ACTION_(NEUROX_ACTION_ZERO_VAR_, DataLoader::InitNetcons);
  NEUROX_REGISTER_ACTION_(NEUROX_ACTION_ZERO_VAR_, DataLoader::Finalize);
  NEUROX_REGISTER_ACTION_(NEUROX_ACTION_MULTIPLE_VARS_, DataLoader::AddSynapse);
  NEUROX_REGISTER_ACTION_(NEUROX_ACTION_MULTIPLE_VARS_, DataLoader::AddNeurons);
  NEUROX_REGISTER_ACTION_(NEUROX_ACTION_MULTIPLE_VARS_,
                          DataLoader::SetMechanisms);
}
