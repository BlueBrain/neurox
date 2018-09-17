#pragma once

#include <vector>
#include "neurox/neurox.h"

namespace neurox {
namespace tools {
namespace linear {

/**
 * @brief The Linear Vector class
 */
template <class Val>
class Vector {
 public:
  Vector() = delete;

  Vector(std::vector<Val*> data, unsigned char* buffer) {
    // needs to be called with placement-new where buffer*
    // is the start of the data structure
    assert(buffer == this);
    n_ = data.size();
    size_t offset = sizeof(Vector<Val>);
    data_ = (Val*)&(buffer[offset]);
    for (int i = 0; i < data.size(); i++)
      memcpy(&data_[i], data[i], sizeof(Val));
  }

  static size_t Size(size_t n) { return sizeof(Vector<Val>) + sizeof(Val) * n; }

  inline Val At(size_t i) {
    assert(i < n_);
    return data_[i];
  }

  inline size_t Count() { return n_; }

 private:
  size_t n_;
  Val* data_;
};  // class Vector
};  // namespace linear
};  // namespace tools
};  // namespace neurox
