#include <numeric>
#include "coreneuron/utils/randoms/nrnran123.h"  //RNG data structures
#include "neurox/neurox.h"

using namespace neurox;
using namespace neurox::tools;

Statistics::SizeInfo::SizeInfo()
    : neuron_id_(0),
      morphologies_(0),
      mechanisms_(0),
      synapses_(0),
      metadata_(0),
      global_vars_(0),
      compartments_count_(0),
      branches_count_(0),
      mechs_instances_count_(0){};

Statistics::SizeInfo::~SizeInfo(){};

double Statistics::SizeInfo::getTotalSize() {
  return morphologies_ + mechanisms_ + synapses_ + metadata_ + global_vars_;
}

Statistics::SizeInfo& Statistics::SizeInfo::operator+=(const SizeInfo& rhs) {
  mechanisms_ += rhs.mechanisms_;
  metadata_ += rhs.metadata_;
  morphologies_ += rhs.morphologies_;
  synapses_ += rhs.synapses_;
  global_vars_ += rhs.global_vars_;
  compartments_count_ += rhs.compartments_count_;
  branches_count_ += rhs.branches_count_;
  mechs_instances_count_ += rhs.mechs_instances_count_;
  return *this;
}

void Statistics::OutputSimulationSize(bool write_to_file) {
  if (hpx_get_my_rank() == 0)
    printf("neurox::tools::Statistics::OutputMechanismsDistribution...\n");

  SizeInfo sim_size;
  sim_size.global_vars_ =
      (double)(sizeof(hpx_t) + sizeof(int) * 2 +
               sizeof(Mechanism) * mechanisms_count_ +
               sizeof(neurox::tools::CmdLineParser) * HPX_LOCALITIES) /
      1024;

  FILE* outstream = stdout;
  if (write_to_file)
    outstream = fopen(string("neurons-memory-consumption.csv").c_str(), "wt");

  fprintf(outstream,
          "gid,compartments,branches,mechs-instances,total-KB,morphologies-KB,"
          "mechanisms-KB,synapses-KB,metadata-KB\n");

  for (int i = 0; i < neurox::neurons_count_; i++) {
    SizeInfo neuron_size;
    hpx_call_sync(neurox::neurons_[i], Statistics::GetNeuronSize, &neuron_size,
                  sizeof(neuron_size));
    fprintf(outstream, "%d,%llu,%llu,%llu,%.1f,%.2f,%.2f,%.2f,%.2f\n",
            neuron_size.neuron_id_, neuron_size.compartments_count_,
            neuron_size.branches_count_, neuron_size.mechs_instances_count_,
            neuron_size.getTotalSize(), neuron_size.morphologies_,
            neuron_size.mechanisms_, neuron_size.synapses_,
            neuron_size.metadata_);
    sim_size += neuron_size;
  }

  printf(
      "- SUM %llu neurons, %llu branches, %llu compartments, %llu mech "
      "instances, %.1f MB\n",
      neurox::neurons_count_, sim_size.branches_count_,
      sim_size.compartments_count_, sim_size.mechs_instances_count_,
      sim_size.getTotalSize() / 1024);
  printf(
      "- AVG per neuron: %.2f branches, %.2f compartments, %.2f mech "
      "instances, %.2f KB\n",
      sim_size.branches_count_ / (double)neurox::neurons_count_,
      sim_size.compartments_count_ / (double)neurox::neurons_count_,
      sim_size.mechs_instances_count_ / (double)neurox::neurons_count_,
      sim_size.getTotalSize() / (double)neurox::neurons_count_);
  printf(
      "- SUM morphologies %.2f MB, mechanisms %.2f MB, synapses %.2f MB, "
      "metadata %.2f MB;\n",
      sim_size.morphologies_ / 1024., sim_size.mechanisms_ / 1024.,
      sim_size.synapses_ / 1024, sim_size.metadata_ / 1024);
  printf(
      "- AVG per neuron: morphologies %.2f KB, mechanisms %.2f KB, synapses "
      "%.2f KB, metadata %.2f KB;\n",
      sim_size.morphologies_ / (double)neurox::neurons_count_,
      sim_size.mechanisms_ / (double)neurox::neurons_count_,
      sim_size.synapses_ / (double)neurox::neurons_count_,
      sim_size.metadata_ / (double)neurox::neurons_count_);
  printf("- Global vars: %.2f KB (Global data %.2f KB * %d localities)\n",
         sim_size.global_vars_, sim_size.global_vars_ / HPX_LOCALITIES,
         HPX_LOCALITIES);
  if (write_to_file) fclose(outstream);
}

hpx_action_t Statistics::GetNeuronSize = 0;
int Statistics::GetNeuronSize_handler() {
  NEUROX_MEM_PIN(Branch);
  assert(local->nt_->end > 0);
  SizeInfo branch_size;
  int n = local->nt_->end;
  if (local->soma_) {
    branch_size.neuron_id_ += local->soma_->gid_;
    branch_size.metadata_ += (double)sizeof(Neuron) / 1024;
    branch_size.synapses_ +=
        (double)(local->soma_->GetSynapsesCount() * sizeof(Neuron::Synapse)) /
        1024;
  }
  branch_size.branches_count_++;
  branch_size.compartments_count_ += n;
  branch_size.morphologies_ +=
      (double)(n * (sizeof(floble_t) * 6)) / 1024;  // a,b,d,v,rhs,area
  branch_size.morphologies_ +=
      local->nt_->_v_parent_index ? (double)(n * sizeof(offset_t)) / 1024 : 0;
  if (local->branch_tree_)
    branch_size.morphologies_ +=
        local->branch_tree_->branches_
            ? (double)(local->branch_tree_->branches_count_ * sizeof(hpx_t)) /
                  1024
            : 0;
  branch_size.metadata_ += (double)sizeof(Branch) / 1024;
  branch_size.metadata_ += (double)sizeof(Memb_list) * mechanisms_count_ / 1024;

  for (int m = 0; m < mechanisms_count_; m++) {
    if (local->mechs_instances_[m].nodecount == 0) continue;

    branch_size.mechs_instances_count_ += local->mechs_instances_[m].nodecount;
    branch_size.mechanisms_ +=
        (double)(sizeof(offset_t) * local->mechs_instances_[m].nodecount) /
        1024;
    if (mechanisms_[m]->data_size_ > 0)
      branch_size.mechanisms_ +=
          (double)(sizeof(floble_t) * mechanisms_[m]->data_size_ *
                   local->mechs_instances_[m].nodecount) /
          1024;
    if (mechanisms_[m]->pdata_size_ > 0)
      branch_size.mechanisms_ +=
          (double)(sizeof(offset_t) * mechanisms_[m]->pdata_size_ *
                   local->mechs_instances_[m].nodecount) /
          1024;
    if (mechanisms_[m]->vdata_size_ > 0)
      branch_size.mechanisms_ +=
          (double)(sizeof(void*) * mechanisms_[m]->vdata_size_ *
                   local->mechs_instances_[m].nodecount) /
          1024;

    int type = mechanisms_[m]->type_;
    if (type == MechanismTypes::kIClamp ||
        type == MechanismTypes::kProbAMPANMDA_EMS ||
        type == MechanismTypes::kProbGABAAB_EMS)
      branch_size.synapses_ += ((double)sizeof(Point_process)) / 1024;
    if (type == MechanismTypes::kStochKv ||
        type == MechanismTypes::kProbAMPANMDA_EMS ||
        type == MechanismTypes::kProbGABAAB_EMS)
      branch_size.synapses_ += ((double)sizeof(nrnran123_State)) / 1024;
  }

  // call the print function in children branches, pass their size to parent
  // branch
  if (local->branch_tree_ && local->branch_tree_->branches_count_ > 0) {
    int branches_count = local->branch_tree_->branches_count_;
    SizeInfo* sub_branches_size = new SizeInfo[branches_count];

    hpx_t* futures = new hpx_t[branches_count];
    void** addrs = new void*[branches_count];
    size_t* sizes = new size_t[branches_count];
    for (offset_t c = 0; c < branches_count; c++) {
      futures[c] = hpx_lco_future_new(sizeof(SizeInfo));
      addrs[c] = &sub_branches_size[c];
      sizes[c] = sizeof(SizeInfo);
      hpx_call(local->branch_tree_->branches_[c], Statistics::GetNeuronSize,
               futures[c]);
    }
    hpx_lco_get_all(branches_count, futures, sizes, addrs, NULL);
    hpx_lco_delete_all(branches_count, futures, NULL);

    delete[] futures;
    delete[] addrs;
    delete[] sizes;

    for (int c = 0; c < branches_count; c++)
      branch_size += sub_branches_size[c];

    delete[] sub_branches_size;
  }

  NEUROX_MEM_UNPIN_CONTINUE(branch_size);  // TODO
}

void Statistics::OutputMechanismsDistribution(bool write_to_file) {
  if (hpx_get_my_rank() == 0)
    printf("neurox::tools::Statistics::OutputMechanismsDistribution...\n");

  unsigned* mechs_count_per_type = new unsigned[mechanisms_count_];
  unsigned* sum_mechs_count_per_type = new unsigned[mechanisms_count_]();
  std::vector<int>* count_per_type = new std::vector<int>[mechanisms_count_];

  unsigned long long total_mechs_instances = 0;
  for (int i = 0; i < neurox::neurons_count_; i++) {
    hpx_call_sync(neurox::neurons_[i],
                  Statistics::GetNeuronMechanismsDistribution,
                  mechs_count_per_type, sizeof(unsigned) * mechanisms_count_);
    for (int m = 0; m < mechanisms_count_; m++) {
      sum_mechs_count_per_type[m] += mechs_count_per_type[m];
      total_mechs_instances += mechs_count_per_type[m];
      count_per_type[m].push_back(mechs_count_per_type[m]);
    }
  }
  delete[] mechs_count_per_type;
  printf("- Total mechs instances: %lld\n", total_mechs_instances);

  FILE* outstream = stdout;
  if (write_to_file)
    outstream = fopen(string("mechs-distribution.csv").c_str(), "wt");
  fprintf(outstream, "mech-type,name,instances,mean-per-neuron,stddev\n");

  for (int m = 0; m < mechanisms_count_; m++) {
    double mean = (double)sum_mechs_count_per_type[m] / neurox::neurons_count_;

    // standard deviation = sqrt ( 1/N * sum (x_i - mean)^2 )
    double std_dev = 0;
    for (int& count : count_per_type[m])
      std_dev += (count - mean) * (count - mean);
    std_dev = sqrt(std_dev / (double)count_per_type[m].size());

    count_per_type[m].clear();

    fprintf(outstream, "%d,%s,%d,%.2f,%.2f\n", mechanisms_[m]->type_,
            mechanisms_[m]->memb_func_.sym, sum_mechs_count_per_type[m], mean,
            std_dev);
  }

  if (write_to_file) fclose(outstream);

  delete[] count_per_type;
  delete[] sum_mechs_count_per_type;
}

hpx_action_t Statistics::GetNeuronMechanismsDistribution = 0;
int Statistics::GetNeuronMechanismsDistribution_handler() {
  NEUROX_MEM_PIN(Branch);
  unsigned mechs_count_per_Type[mechanisms_count_];
  for (int m = 0; m < mechanisms_count_; m++)
    mechs_count_per_Type[m] = local->mechs_instances_[m].nodecount;

  // call the function on children branches, pass their size to parent branch
  if (local->branch_tree_ && local->branch_tree_->branches_count_ > 0) {
    int branches_count = local->branch_tree_->branches_count_;

    unsigned** mechs_count_per_type_child = new unsigned*[branches_count];
    for (int c = 0; c < branches_count; c++)
      mechs_count_per_type_child[c] = new unsigned[mechanisms_count_];

    hpx_t* futures = new hpx_t[branches_count];
    void** addrs = new void*[branches_count];
    size_t* sizes = new size_t[branches_count];
    for (offset_t c = 0; c < branches_count; c++) {
      futures[c] = hpx_lco_future_new(sizeof(unsigned[mechanisms_count_]));
      addrs[c] = mechs_count_per_type_child[c];
      sizes[c] = sizeof(unsigned[mechanisms_count_]);
      hpx_call(local->branch_tree_->branches_[c],
               Statistics::GetNeuronMechanismsDistribution, futures[c]);
    }
    hpx_lco_get_all(branches_count, futures, sizes, addrs, NULL);
    hpx_lco_delete_all(branches_count, futures, NULL);

    delete[] futures;
    delete[] addrs;
    delete[] sizes;

    for (int c = 0; c < branches_count; c++)
      for (int m = 0; m < mechanisms_count_; m++)
        mechs_count_per_Type[m] += mechs_count_per_type_child[c][m];

    for (int c = 0; c < branches_count; c++)
      delete[] mechs_count_per_type_child[c];
    delete[] mechs_count_per_type_child;
  }
  NEUROX_MEM_UNPIN_CONTINUE(mechs_count_per_Type);  // TODO
}

void Statistics::RegisterHpxActions() {
  wrappers::RegisterZeroVarAction(GetNeuronSize, GetNeuronSize_handler);
  wrappers::RegisterZeroVarAction(
      Statistics::GetNeuronMechanismsDistribution,
      Statistics::GetNeuronMechanismsDistribution_handler);
}
