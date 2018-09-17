#pragma once

#include "neurox/neurox.h"

namespace neurox {
namespace tools {
namespace linear {

/**
 * @brief The Linear Vector class
 */
template <class T>
class Vector {
 public:
  // needs to be called with placement-new where buffer*
  // is the start of the data structure
  Vector(size_t n, T* data, char* buffer) : n_(n) {
    assert(buffer == this);
    size_t offset = sizeof(n) + sizeof(data);
    data_ = (T*)&(buffer[offset]);
    memcpy(data_, data, sizeof(T) * n);
  }

  size_t Size(size_t n) { return sizeof(Vector<T>) + sizeof(T) * n; }

  inline T At(size_t i) {
    assert(i < n_);
    return data_[i];
  }

  inline size_t Count() { return n_; }

 private:
  size_t n_;
  T* data_;
};  // class Vector
};  // namespace linear
};  // namespace tools
};  // namespace neurox
