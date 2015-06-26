#pragma once

#include "neurox/neurox.h"

#include <vector>

namespace neurox {

class NetConX;

/// hard-coded mechanism types
enum MechanismTypes {
  kIClamp = 7,
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

  Mechanism(const int type, const short dataSize, const short pdataSize,
            const char isArtificial, char pntMap, const char isIon,
            const short int symLengh, const char *sym, Memb_func &memb_func,
            const short int dependenciesCount = 0,
            const int *dependencies = nullptr,
            const short int successorsCount = 0,
            const int *successors = nullptr);

  int type;
  short dataSize, pdataSize, vdataSize;
  short successorsCount;    ///> number of mechanisms succedding this one on a
                            /// parallel execution
  short dependenciesCount;  ///> number of mechanisms it depends on
  short symLength;          ///> length of the name of the mechanism;
  char pntMap, isArtificial;
  char isIon;
  int *dependencies;  ///> mechanism id for dependency mechanisms
  int *successors;    ///> mechanism id for successors mechanisms

  int dependencyIonIndex;  ///> index of parent ion (if any)

  // from memb_func.h (before after functions not used on BBP models)
  Memb_func membFunc;
  mod_f_t BAfunctions[BEFORE_AFTER_SIZE];  ///>mechanism functions
  pnt_receive2_t pnt_receive;
  pnt_receive2_t pnt_receive_init;
  bbcore_read_t nrn_bbcore_read;

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
    kNetReceiveInit = 18
  };

  Mechanism::IonTypes GetIonIndex();

  void CallModFunction(const void *branch,
                       const Mechanism::ModFunctions functionId,
                       const NetConX *netcon = NULL,  // for net_receive only
                       const floble_t tt = 0);        // for net_receive only
 private:
  void RegisterBeforeAfterFunctions();  ///> register Before-After functions
};
};
