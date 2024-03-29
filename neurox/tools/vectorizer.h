/*
# =============================================================================
# Copyright (c) 2015 - 2021-2022 Blue Brain Project/EPFL
#
# See top-level LICENSE file for details.
# =============================================================================.
*/
#pragma once

#include "neurox/neurox.h"

#include <deque>
#include <map>
#include <memory>
#include <set>
#include <vector>

using namespace std;

namespace neurox {
namespace tools {

/**
 * @brief The Vectorizer class
 * converts code from AoS to SoA structure
 */
class Vectorizer {
 public:
  Vectorizer() = delete;
  ~Vectorizer() = delete;

  /// memory alignment used for vectorized data structures
  static const int kMemoryAlignment = 2 * sizeof(double);

  /// padding added to vectorized on memory-aligned data-structures
  static const int kSOAPadding = 4;

  /// converts a branch from AoS to SoA
  static void ConvertToSOA(Branch *b);

  /// get size of datastructure with alignment-based padding
  static size_t SizeOf(size_t size) {
    return coreneuron::soa_padded_size<kSOAPadding>(size, LAYOUT);
  }

  /// Call function 'f' on a vectorized way
  static void CallVecFunction(cvode_f_t, NrnThread *, Memb_list *, int);

  /// returns a copy of Memb_list of a branch, sorted by no-cap and cap
  static void GroupBranchInstancesByCapacitors(
      const Branch *branch,                        // in
      Memb_list **ml_no_capacitors_ptr = nullptr,  // out (optional)
      Memb_list **ml_capacitors_ptr = nullptr,     // out (optional)
      std::set<int> *capacitors_ids_ptr = nullptr  // in (optional)
  );

  /// creates branch->
  static void CreateMechInstancesThreads(Branch *branch);

  // C++11 does not support memory-aligned new[]/delete, this is a work around

  /// memory-aligned memory allocation
  template <typename T>
  static T *New(size_t count) {
    return count == 0 ? nullptr
                      : (T *)coreneuron::ecalloc_align(
                            SizeOf(count), kMemoryAlignment, sizeof(T));
  }

  template <typename T>
  static T *New(size_t count, void *buffer, size_t buffer_size,
                size_t &offset) {
    if (count == 0) return nullptr;

    // insteaf of allocating memory, send next position in buffer available
    if (buffer) {
      assert(offset < buffer_size);
      void *pos = &((unsigned char *)buffer)[offset];
      offset += SizeOf(sizeof(T) * count);
      return (T *)pos;
    }
    return New<T>(count);
  }

  /// delete for the New method
  template <typename T>
  static void Delete(T *ptr) {
    free(ptr);
    ptr = nullptr;  // delete[] (ptr);
  }

 private:
};

};  // namespace tools
};  // namespace neurox
