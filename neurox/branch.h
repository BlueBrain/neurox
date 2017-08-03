#pragma once

#include "neurox/neurox.h"

#include <functional>  //std::greater_equal
#include <map>
#include <new>  //placement new
#include <queue>
#include <vector>

using namespace neurox;

namespace neurox {

class Neuron;

/**
 * @brief The Branch class
 * Represents a branch as a continuous set of compartments (nodes)
 * and a bifurcation made of branches branchs spreading out from the last node;
 * Also represents mechanisms as a list of mechanisms and their applications
 * to several nodes of this branch.
 */
class Branch {
 public:
  Branch() = delete;  ///> no constructor, build using hpx init function instead
  static void* operator new(size_t bytes, void* addr);
  static void operator delete(void* worker);

  Branch(offset_t n, int nrn_thread_id, int threshold_v_offset,
         hpx_t branch_hpx_addr, floble_t* data, size_t data_count,
         offset_t* pdata, size_t pdata_count, offset_t* instances_count,
         size_t recv_mechs_count, offset_t* nodes_indices,
         size_t nodes_indices_count, hpx_t top_branch_addr, hpx_t* branches,
         size_t branches_count, offset_t* p, size_t p_count,
         floble_t* vecplay_t, size_t vecplay_t_Count, floble_t* vecplay_y,
         size_t vecplay_y_count, PointProcInfo* vecplay_ppi,
         size_t vecplay_ppi_count, NetconX* netcons_, size_t netcons_count,
         neuron_id_t* netcons_pre_ids, size_t netcons_pre_ids_count,
         floble_t* weights, size_t weights_count,
         unsigned char* vdata_serialized, size_t vdata_serialized_count);
  ~Branch();

  NrnThread* nt_;               ///> compartments metadata
  Memb_list* mechs_instances_;  ///> Arrays of mechanism instances
  Neuron* soma_;         ///> if top branch, it's populated, otherwise NULL
  floble_t* thvar_ptr_;  ///> pointer to var holding AP threshold var,if any

  class MechanismsGraph {
   public:
    MechanismsGraph();
    ~MechanismsGraph();

    /// Launches HPX-threads for dorment mechs-graph
    void InitMechsGraph(hpx_t branch_hpx_addr);

    hpx_t* mechs_lcos_;  ///> HPX address of the and-gate of each mechanism
    hpx_t end_lco_;      ///> represents the bottom of the graph
    hpx_t graph_lco_;  ///> controls all active threads on the mechanisms graph

    /// init function for hpx_reduce of mechanisms graphs
    static hpx_action_t Init;

    /// reduce opreation for hpx_reduce of mechanisms graphs
    static hpx_action_t Reduce;

    /// represents the action of the nodes in the mechanisms graph
    static hpx_action_t MechFunction;

    static int MechFunction_handler(const int* mech_type_ptr, const size_t);
    static void Init_handler(Mechanism::ModFunctions* func_ptr, const size_t);
    static void Reduce_handler(Mechanism::ModFunctions* lhs,
                               const Mechanism::ModFunctions* rhs,
                               const size_t);

    // for current function accumulation of shadow arrays
    hpx_t rhs_d_mutex_;
    hpx_t i_didv_mutex_[Mechanism::IonTypes::kSizeWriteableIons];
    static void AccumulateRHSandD(NrnThread* nt, Memb_list* ml, int,
                                  void* args);
    static void AccumulateIandDIDV(NrnThread* nt, Memb_list* ml, int,
                                   void* args);

  } * mechs_graph_;  ///> parallel computation graph of mechanisms

  class BranchTree {
   public:
    BranchTree() = delete;
    BranchTree(hpx_t top_branch_addr, hpx_t* branches, size_t branches_count);
    ~BranchTree();

    hpx_t top_branch_addr_;  ///> hpx address of the some branch
    hpx_t* branches_;        ///> hpx address of children branches
    size_t branches_count_;  ///> number of branches (>0)

    ///  size of futures arrays (used in Gaussian elimination and AP threshold
    static constexpr size_t kFuturesSize = 7;

    /// LCO to to communicate variables with parent
    hpx_t with_parent_lco_[kFuturesSize];

    /// LCO to communicate variables with children (NULL if no children)
    hpx_t (*with_children_lcos_)[kFuturesSize];

    static hpx_action_t InitLCOs;  ///> Initializes neuronTree
    static int InitLCOs_handler();
  } * branch_tree_;  ///> represents the tree structure (or NULL if none)

  /// map of incoming netcons per pre-synaptic gid
  std::map<neuron_id_t, std::vector<NetconX*> > netcons_;

  /// priority queue of incoming events sorted per delivery time
  std::priority_queue<TimedEvent, std::vector<TimedEvent>,
                      std::greater_equal<TimedEvent> >
      events_queue_;

  /// mutex to protect the memory access to eventsQueue
  hpx_t events_queue_mutex_;

  static hpx_action_t Init;  ///> Initializes the diagonal matrix and branching

  static hpx_action_t InitSoma;  ///> Initializes soma information
  static hpx_action_t Clear;  ///> deletes all data in branch and sub-branches
  static hpx_action_t AddSpikeEvent;  ///>add incoming synapse to queue
  /// update maximum time allowed based on received dependency info
  static hpx_action_t UpdateTimeDependency;
  static hpx_action_t Finitialize;  ///> finitialize.c::finitialize()
  static hpx_action_t BackwardEuler;
  static hpx_action_t BackwardEulerOnLocality;
  static hpx_action_t ThreadTableCheck;

  void CallModFunction(const Mechanism::ModFunctions functionId);
  void
  InitVecPlayContinous();  ///> start NetEvents and PlayVect on events queue
  void AddEventToQueue(floble_t t, Event* e);
  void DeliverEvents(floble_t t);
  void FixedPlayContinuous();
  void SetupTreeMatrix();
  void SolveTreeMatrix();
  void Finitialize2();
  void BackwardEulerStep();

  static void RegisterHpxActions();  ///> Register all HPX actions

 private:
  static int Init_handler(const int, const void* [], const size_t[]);
  static int InitSoma_handler(const int, const void* [], const size_t[]);
  static int Clear_handler();
  static int AddSpikeEvent_handler(const int, const void* [], const size_t[]);
  static int UpdateTimeDependency_handler(const int, const void* [],
                                          const size_t[]);
  static int Finitialize_handler();
  static int BackwardEuler_handler(const int*, const size_t);
  static int BackwardEulerOnLocality_handler(const int*, const size_t);
  static int ThreadTableCheck_handler();
};

};  // namespace