#pragma once

#include "neurox/neurox.h"

#include <deque>
#include <map>
#include <memory>
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
  static void ConvertToSOA(Branch* b);

  /// get size of datastructure with alignment-based padding
  static size_t SizeOf(size_t size) {
    return coreneuron::soa_padded_size<kSOAPadding>(size, LAYOUT);
  }

  // C++11 does not support memory-aligned new[]/delete, this is a work around
  template <typename T>
  static T* New(size_t count) {
    return (T*)coreneuron::ecalloc_align(SizeOf(count), kMemoryAlignment,
                                         sizeof(T));
  }

  template <typename T>
  static void Delete(T* ptr) {
    free(ptr);
    ptr = nullptr;  // delete[] (ptr);
  }

 private:
};

};  // Tools
};  // NeuroX
