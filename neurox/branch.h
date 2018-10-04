#pragma once

#include "neurox/neurox.h"

#include <functional>  //std::greater_equal
#include <map>
#include <new>  //placement new
#include <queue>
#include <vector>

using namespace neurox;
using namespace tools;

namespace neurox {

// Fwd declarations
namespace interpolators {
class Interpolator;
}

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

  Branch(offset_t n, int nrn_thread_id, int threshold_v_offset, floble_t* data,
         size_t data_count, offset_t* pdata, size_t pdata_count,
         offset_t* instances_count, size_t recv_mechs_count,
         offset_t* nodes_indices, size_t nodes_indices_count,
         hpx_t top_branch_addr, hpx_t* branches, size_t branches_count,
         offset_t* p, size_t p_count, floble_t* vecplay_t,
         size_t vecplay_t_Count, floble_t* vecplay_y, size_t vecplay_y_count,
         PointProcInfo* vecplay_ppi, size_t vecplay_ppi_count,
         NetconX* netcons_, size_t netcons_count, neuron_id_t* netcons_pre_ids,
         size_t netcons_pre_ids_count, floble_t* weights, size_t weights_count,
         unsigned char* vdata_serialized, size_t vdata_serialized_count,
         neuron_id_t soma_gid, floble_t soma_ap_threshold);
  ~Branch();

  /// if using cache-effificent representation, all data is serialized here
  unsigned char* buffer_;
  size_t buffer_size_;

  /// clears a given Memb_list structure
  static void ClearMembList(Memb_list*&);

  /// clears a given NrnThread structure
  static void ClearNrnThread(NrnThread*&);

  /// Compartments, Weights and other Coreneuron data
  NrnThread* nt_;

  /// Array of mechanisms instances
  Memb_list* mechs_instances_;

  /** Array of arguments for parallel-threaded execution of
   * mechanisms instances, or nullptr when not applicable */
  Mechanism::MembListThreadArgs* mechs_instances_parallel_;

  /// if top branch, refers to soma/neuron info, otherwise is NULL
  Neuron* soma_;

  /// pointer to var holding AP threshold var,if any
  floble_t* thvar_ptr_;

  //// graph-based execution of mechanisms
  class MechanismsGraph {
   public:
    MechanismsGraph();
    ~MechanismsGraph();

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

  }* mechs_graph_;  ///> parallel computation graph of mechanisms

  class BranchTree {
   public:
    BranchTree() = delete;
    BranchTree(hpx_t top_branch_addr, hpx_t* branches, size_t branches_count);
    ~BranchTree();

    hpx_t top_branch_addr_;  ///> hpx address of the some branch
    hpx_t* branches_;        ///> hpx address of children branches
    size_t branches_count_;  ///> number of branches (>0)

    /// size of futures arrays (used in Gaussian elimination and AP threshold
    /// channel 0: previous-step RHS from children (after ForwardSubstitution)
    /// channel 1: D+RHS from children (BackwardTriangulation)
    /// channel 2: RHS from parent (after ForwardSubstitution)
    static constexpr size_t kFuturesSize = 3;

    /// value of A[0] of all children
    floble_t* children_a_;

    /// value of B[0] of all children
    floble_t* children_b_;

    /// value of V[0] of all children
    floble_t* children_v_;

    /// value of RHS[0] of all children
    floble_t* children_rhs_;

    /// value of V[n-1] from parent
    floble_t parent_v_;

    /// value of RHS[n-1] from parent
    floble_t parent_rhs_;

    /// LCO to to communicate variables with parent
    hpx_t with_parent_lco_[kFuturesSize];

    /// LCO to communicate variables with children (NULL if no children)
    hpx_t (*with_children_lcos_)[kFuturesSize];

    static hpx_action_t InitLCOs;  ///> Initializes neuronTree
    static int InitLCOs_handler();
  }* branch_tree_;  ///> represents the tree structure (or NULL if none)

  /// map of incoming netcons per pre-synaptic gid
  std::map<neuron_id_t, std::vector<NetconX*> > netcons_;

  linear::Map<neuron_id_t, NetconX>* netcons_linear_;

  /// priority queue of incoming events sorted per delivery time
  std::priority_queue<TimedEvent, std::vector<TimedEvent>,
                      std::greater_equal<TimedEvent> > events_queue_;

  /// for a neuron_id_t, gives the next Event* on the list
  /// (we use Event* cause it points to data structs thar are linearized)
  linear::PriorityQueue<neuron_id_t, TimedEvent>* events_queue_linear_;

  /// mutex to protect the memory access to eventsQueue
  hpx_t events_queue_mutex_;

  static hpx_action_t Init;  ///> Initializes the diagonal matrix and branching

  static hpx_action_t InitMechanismsGraph;  ///> Initializes mechanisms graph
  static hpx_action_t InitMechParallelism;  ///> Initializes M.I.P.
  static hpx_action_t Initialize;  ///> Initializes interpolator for this neuron
  static hpx_action_t Clear;  ///> deletes all data in branch and sub-branches
  static hpx_action_t AddSpikeEvent;          ///>add incoming synapse to queue
  static hpx_action_t AddSpikeEventLocality;  ///>add incoming synapse to queue
  static hpx_action_t SetSyncStepTrigger;

  /// update maximum time allowed based on received dependency info
  static hpx_action_t ThreadTableCheck;

  void CallModFunction(const Mechanism::ModFunctions functionId,
                       Memb_list* other_ml = nullptr);

  /// start NetEvents and PlayVect on events queue
  void InitVecPlayContinous();

  static void CoreneuronNetSend(void**, int, NrnThread*, int, int, double,
                                double);
  static void CoreneuronNetEvent(NrnThread*, int, int, double);

  void AddEventToQueue(floble_t t, Event* e);
  void DeliverEvents(floble_t t);
  void FixedPlayContinuous(double);
  void FixedPlayContinuous();
  void SetupTreeMatrix();

  static void RegisterHpxActions();  ///> Register all HPX actions

  /// Interpolator for fixed- or variable- stepping
  interpolators::Interpolator* interpolator_;

 private:
  static int Init_handler(const int, const void * [], const size_t[]);
  static int InitMechanismsGraph_handler();
  static int InitMechParallelism_handler();
  static int Initialize_handler();
  static int Clear_handler();
  static int AddSpikeEvent_handler(const int, const void * [], const size_t[]);
  static int AddSpikeEventLocality_handler(const int, const void * [],
                                           const size_t[]);
  static int ThreadTableCheck_handler();
  static int SetSyncStepTrigger_handler(const hpx_t*, const size_t);
};

};  // namespace neurox
