#pragma once

#include "neurox/neurox.h"

#include <deque>
#include <map>
#include <memory>
#include <set>
#include <vector>

#define NEUROX_INPUT_DATALOADER_OUTPUT_EXTERNAL_NETCONS true
#define NEUROX_INPUT_DATALOADER_OUTPUT_CORENEURON_COMPARTMENTS false

using namespace std;

namespace neurox {

namespace input {

/**
 * @brief The NrxSetup class
 * Converts CoreNeuron data structures to HPX data structures
 */
class DataLoader {
 public:
  DataLoader() = delete;
  ~DataLoader() = delete;

  /// Copies Coreneuron data structs to HPX memory
  static void loadData(int argc, char **argv);

  /// Calls coreneuron nrn_init_and_load_data
  static void LoadCoreneuronData(int argc, char **argv,
                                 bool nrnmpi_under_nrncontrol = false,
                                 bool run_setup_cleanup = false);

  /// removes all data structures loaded for coreneuron
  static void CleanCoreneuronData(const bool clean_ion_global_map = true);

  static void RegisterHpxActions();  ///> Register HPX actions

  static hpx_action_t Init;
  static hpx_action_t InitMechanisms;
  static hpx_action_t InitNeurons;
  static hpx_action_t SetNeurons;
  static hpx_action_t InitNetcons;
  static hpx_action_t Finalize;
  static hpx_action_t FilterRepeatedAndLinearizeContainers;

  /// holds information about offsets of instances of mechanisms
  struct IonInstancesInfo {
   public:
    int mech_type;
    offset_t data_start;   ///> beginning offset of instances in data
    offset_t data_end;     ///> end offset of instances in data
    vector<int> node_ids;  ///> compartments ids for each instance
  };

  /// mutex controlling multi-threaded write of data structs
  static hpx_t locality_mutex_;

  static int HardCodedVdataSize(int type);
  static int HardCodedVdataCount(int type, char pnt_map);
  static int HardCodedPntProcOffsetInPdata(int type);
  static int HardCodedPntProcOffsetInVdata(int type);
  static int HardCodedQueueItemOffsetInPdata(int type);
  static int HardCodedPPtrOffsetInPdata(int type);
  static int HardCodedRNGOffsetInPdata(int type);
  static int HardCodedRNGOffsetInVdata(int type);
  static bool HardCodedMechanismHasNoInstances(int index);
  static bool HardCodedMechanismForCoreneuronOnly(int index);
  static bool HardCodedEventIsDiscontinuity(Event *);

 private:
  /// pointer of netcons.dot file
  static FILE *file_netcons_;

  /// temporary hpx address of neurons read by this locality
  static std::vector<hpx_t> *my_neurons_addr_;

  /// temporary gid of neurons read by this locality
  static std::vector<int> *my_neurons_gids_;

  /// gids of all neurons in the system
  static std::vector<int> *all_neurons_gids_;

  /// pointer to load balancing instantiated class (if any)
  static tools::LoadBalancing *load_balancing_;

  /// compute execution time of a given subsection
  static double BenchmarkSubSection(int N,
                                    const deque<Compartment *> &sub_section,
                                    vector<DataLoader::IonInstancesInfo> &);

  static int GetNumberOfInstanceCompartments(
      const deque<Compartment *> &compartments);

  static hpx_t CreateBranch(
      const int nrn_thread_id, hpx_t soma_branch_addr,
      const deque<Compartment *> &all_compartments,
      Compartment *top_compartment,
      vector<DataLoader::IonInstancesInfo> &ions_instances_info,
      double neuron_time, int thvar_index = -1 /*AIS*/,
      floble_t ap_threshold = 0 /*AIS*/, int assigned_locality = -1);

  static neuron_id_t GetNeuronIdFromNrnThreadId(int nrn_id);
  static void getMechTypeAndInstanceForBranch(int &mech_type,
                                              int &mech_instance);

  static int GetBranchData(
      const deque<Compartment *> &compartments, vector<floble_t> &data,
      vector<offset_t> &pdata, vector<unsigned char> &vdata,
      vector<offset_t> &p, vector<offset_t> &instances_count,
      vector<offset_t> &nodes_indices, int N,
      vector<DataLoader::IonInstancesInfo> &ions_instances_info,
      vector<map<int, int>> *mech_instance_map = NULL);

  static void GetVecPlayBranchData(
      const deque<Compartment *> &compartments,
      vector<floble_t> &vecplay_t_data, vector<floble_t> &vecplay_y_data,
      vector<PointProcInfo> &vecplay_info,
      vector<map<int, int>> *mech_instance_map = NULL);

  static void GetNetConsBranchData(
      const deque<Compartment *> &compartments, vector<NetconX> &branch_netcons,
      vector<neuron_id_t> &branch_netcons_pre_id,
      vector<floble_t> &branch_netcons_args,
      vector<map<int, int>> *mech_instance_map = NULL);

  static void GetSubSectionFromCompartment(deque<Compartment *> &sub_section,
                                           Compartment *top_compartment);

  static void GetMechInstanceMap(const deque<Compartment *> &compartments,
                                 vector<map<int, int>> &mechs_instance_map);

  static void SetMechanisms2(const int mechs_count, const int *mechs_ids,
                             const int *dependencies_count,
                             const int *dependencies,
                             const int *successors_count,
                             const int *successors);  ///> Set Mechanisms

  static void PrintSubClustersToFile(FILE *file_compartments,
                                     Compartment *top_compartment);

  static PointProcInfo GetPointProcInfoFromDataPointer(NrnThread *nt,
                                                       double *pd, size_t size);

  static hpx_action_t AddSynapse;
  static hpx_action_t AddNeurons;
  static hpx_action_t SetMechanisms;

  static int CreateNeuron(int neuron_idx, void *);
  static int GetMyNrnThreadsCount();
  static int AddSynapse_handler(const int, const void *[], const size_t[]);
  static int FilterRepeatedAndLinearizeContainers_handler();
  static int AddNeurons_handler(const int, const void *[], const size_t[]);
  static int SetMechanisms_handler(const int, const void *[], const size_t[]);
  static int Init_handler();
  static int InitMechanisms_handler();
  static int InitNeurons_handler();
  static int SetNeurons_handler();
  static int InitNetcons_handler();
  static int Finalize_handler();
};

};  // namespace input
};  // namespace neurox
