#pragma once

#include "neurox/neurox.h"

#include <vector>

using namespace coreneuron;

namespace neurox {

class NetconX;

/// TODO: hard-coded mechanism types
enum MechanismTypes {
  kCapacitance = 3,
  kIClamp = 7,
  kExpSyn = 9,
  kCaDynamics_E2 = 36,
  kProbAMPANMDA_EMS = 156,
  kProbGABAAB_EMS = 158,
  kStochKv = 176,
  kStochKv3 = 175,
  kBinReportHelper = 26,
  kBinReports = 27,
  kCoreConfig = 56,
  kMemUsage = 131,
  kProfileHelper = 162,
  kVecStim = 181,
  kNetStim = 18,
  kGluSynapse = 63,
  kALU = 24,
  kInhPoissonStim = 152,
  kPatternStim = 23
};

/**
 * @brief The Mechanisms class
 * Stores the unique metadata of each mechanism
 */
class Mechanism {
 public:
  /// unique id identifying ion mechanisms
  /// TODO these types could be merged with the above?!
  enum IonTypes {
    kNa = 0,
    kK = 1,
    kCa = 2,
    kSizeWriteableIons = 3,
    kTTX = 3,
    kSizeAllIons = 4,
    kNoIon = 9
  };

  Mechanism() = delete;
  ~Mechanism();

  Mechanism(const int type_, const short int data_size,
            const short int pdata_size, const char is_artificial, char pnt_map,
            const char is_ion, const short int sym_length, const char *sym,
            coreneuron::Memb_func &memb_func_,
            const short int dependencies_count = 0,
            const int *dependencies_ = nullptr,
            const short int successors_count = 0,
            const int *successors_ = nullptr);

  int type_;
  short data_size_, pdata_size_, vdata_size_;
  short successors_count_;    ///> number of mechanisms succedding this one on a
                              /// parallel execution
  short dependencies_count_;  ///> number of mechanisms it depends on
  short sym_length_;          ///> length of the name of the mechanism;
  char pnt_map_, is_artificial_;
  char is_ion_;
  int *dependencies_;  ///> mechanism id for dependency mechanisms
  int *successors_;    ///> mechanism id for successors mechanisms

  int dependency_ion_index_;  ///> index of parent ion (if any)

  // from memb_func.h (before after functions not used on BBP models)
  Memb_func memb_func_;
  mod_f_t before_after_functions_[BEFORE_AFTER_SIZE];  ///>mechanism functions
  pnt_receive_t pnt_receive_;
  bbcore_read_t nrn_bbcore_read_;

  // CVODES-specific
  cvode_f_t ode_spec_;
  cvode_f_t ode_matsol_;
  mod_f_t div_capacity_;
  mod_f_t mul_capacity_;

  /// State variables info (used by CVODES only)
  class StateVars {
   public:
    StateVars();
    StateVars(int count, int *offsets, int *dv_offsets);
    ~StateVars();
    int count_;         ///> number of cvode state variables
    int *var_offsets_;  ///>offset of state vars in ml->data
    int *dv_offsets_;   ///> offset of dx/dV for state vars
  } * state_vars_;

  enum ModFunctions {
    // BA functions start here (of size BEFORE_AFTER_SIZE)
    kBeforeInitialize = 0,
    kAfterInitialize = 1,
    kBeforeBreakpoint = 2,
    kAfterSolve = 3,
    kBeforeStep = 4,
    // memb_func functions start here
    kAlloc = 5,
    kCurrent = 6,
    kState = 7,
    kJacob = 8,
    kInitialize = 9,
    kDestructor = 10,
    kThreadMemInit = 11,
    kThreadCleanup = 12,
    kThreadTableCheck = 13,
    kSetData = 14,
    // capacitance functions start here
    kCurrentCapacitance = 15,  // not in mod files, it's in capac.c
    kJacobCapacitance = 16,
    // net_receive
    kNetReceive = 17,
    kNetReceiveInit = 18,
    // CVODE-specific methods
    kODESpec = 19,
    kODEMatsol = 20,
    kDivCapacity = 21,
    kMulCapacity = 22
  };

  Mechanism::IonTypes GetIonIndex();

  void CallModFunction(
      const void *branch, const Mechanism::ModFunctions function_id,
      Memb_list *other_ml = nullptr,    // user-provided Memb_list (if any)
      const NetconX *netcon = nullptr,  // net_receive only
      const floble_t tt = 0);           // net_receive only

  /// argument to mech-instances threaded calls
  typedef struct MembListThreadArgsStruct {
    Memb_func *memb_func;   ///> Pointer to Memb_func in Branch
    NrnThread *nt;          ///> Pointer to NrnThread in Branch
    Memb_list *ml_state;    ///> Array of MembList per thread for state func
    Memb_list *ml_current;  ///> Array of Memblist per phread for current func
    int ml_state_count;     //>  size of array *ml_state
    int ml_current_count;   ///> size of array *ml_current
    int mech_type;          ///> mechanism type

    // parallel execution of current only:
    bool requires_shadow_vectors;
    mod_acc_f_t acc_rhs_d;
    mod_acc_f_t acc_di_dv;
    void *acc_args;
  } MembListThreadArgs;

  /// counter for execution time allocated to mechanisms execution
  static double time_spent_in_mechs_;
  static hpx_t time_spent_in_mechs_mutex_;

 private:
  /// thread function for an instance-parallel call of current function
  static inline int ModFunctionCurrentThread(int i, void *arg_ptr);

  /// thread function for an instance-parallel call of state function
  static inline int ModFunctionStateThread(int i, void *arg_ptr);

  /// register Before-After functions
  void RegisterBeforeAfterFunctions();
};
};  // namespace neurox
