#pragma once

#include <vector>
#include "neurox/neurox.h"

namespace neurox {

/**
 * @brief The Mechanisms class
 * Stores the unique metadata of each mechanism
 */
class Mechanism {
 public:
  enum Ion {
    na = 0,
    k = 1,
    ca = 2,
    size_writeable_ions = 3,
    ttx = 3,
    size_all_ions = 4,
    no_ion = 9
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
                            ///parallel execution
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

  enum ModFunction {
    // BA functions start here (of size BEFORE_AFTER_SIZE)
    before_initialize = 0,
    after_initialize = 1,
    before_breakpoint = 2,
    after_solve = 3,
    before_step = 4,
    // memb_func functions start here
    alloc = 5,
    current = 6,
    state = 7,
    jacob = 8,
    initialize = 9,
    destructor = 10,
    threadMemInit = 11,
    threadCleanup = 12,
    threadTableCheck = 13,
    setData = 14,
    // capacitance functions start here
    currentCapacitance = 15,  // not in mod files, it's in capac.c
    jacobCapacitance = 16,
    // net_receive
    netReceive = 17,
    netReceiveInit = 18
  };

  Mechanism::Ion GetIonIndex();

  void CallModFunction(const void *branch,
                       const Mechanism::ModFunction functionId,
                       const NetConX *netcon = NULL,  // for net_receive only
                       const floble_t tt = 0);        // for net_receive only
 private:
  void RegisterBeforeAfterFunctions();  ///> register Before-After functions
};
};
