#pragma once

#include "neurox/neurox.h"

#include <vector>

namespace neurox {

class NetconX;

/// hard-coded mechanism types
enum MechanismTypes {
  kCapacitance = 3,
  kIClamp = 7,
  kExpSyn = 9,
  kProbAMPANMDA_EMS = 137,
  kProbGABAAB_EMS = 139,
  kStochKv = 151
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
            Memb_func &memb_func_, const short int dependencies_count = 0,
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
  pnt_receive2_t pnt_receive_;
  pnt_receive2_t pnt_receive_init_;
  bbcore_read_t nrn_bbcore_read_;

  //CVODES-specific
  ode_spec1_f_t ode_spec_;
  ode_matsol1_f_t ode_matsol_;

  /// State variables info (used by CVODES only)
  class StateVars {
   public:
    StateVars();
    StateVars(short count, short *offsets, short *dv_offsets);
    ~StateVars();
    short count_;        ///> number of cvode state variables
    short *var_offsets_;     ///>offset of state vars in ml->data
    short *dv_offsets_;  ///> offset of dx/dV for state vars
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
    kODEMatsol = 20
  };

  Mechanism::IonTypes GetIonIndex();

  void CallModFunction(const void *branch,
                       const Mechanism::ModFunctions function_id,
                       const NetconX *netcon = NULL,  // for net_receive only
                       const floble_t tt = 0);        // for net_receive only
 private:
  void RegisterBeforeAfterFunctions();  ///> register Before-After functions
};
};
