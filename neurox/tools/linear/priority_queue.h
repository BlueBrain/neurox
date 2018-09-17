#pragma once

#include "neurox/neurox.h"

namespace neurox {
namespace tools {
namespace linear {

/**
 * @brief The Linear (max fixed size) PriorityQueue class
 */
template <class Key, class Val>
class PriorityQueue {
 public:
  // needs to be called with placement-new where buffer*
  // is the start of the data structure
  PriorityQueue() = delete;

  PriorityQueue(size_t keys_count, Key* keys, size_t* max_vals_per_key,
                unsigned char* buffer)
      : keys_count_(keys_count) {
    assert(buffer == this);

    size_t offset = sizeof(size_t);

    // keys_
    keys_ = (Key*)&(buffer[offset]);
    memcpy(keys_, keys, sizeof(Key) * keys_count);
    offset += sizeof(Key) * keys_count;

    // count of values per key (size of circular array)
    vals_per_key_ = (size_t*)&(buffer[offset]);
    memcpy(vals_per_key_, max_vals_per_key, sizeof(size_t) * keys_count);
    offset += sizeof(size_t) * keys_count;

    // current offset of each circular array
    offsets_per_key_ = (size_t*)&(buffer[offset]);
    for (int i = 0; i < keys_count; i++) offsets_per_key_[i] = 0;
    offset += sizeof(size_t) * keys_count;

    // values array of pointers
    vals_ = (Val**)&(buffer[offset]);
    offset += sizeof(Val*) * keys_count;

    // values
    for (int i = 0; i < keys_count; i++) {
      vals_[i] = (Val*)&(buffer[offset]);
      Val dummy_val;
      for (int j = 0; j < max_vals_per_key[i]; j++)
        memcpy(&(vals_[i][j]), &dummy_val, sizeof(Val));
      i++;
    }
  }

  static size_t Size(size_t keys_count, size_t* max_vals_per_key) {
    size_t size = sizeof(PriorityQueue<Key, Val>);
    size += sizeof(size_t) * keys_count;  // keys
    size += sizeof(size_t) * keys_count;  // vals per key
    size += sizeof(size_t) * keys_count;  // offsets per key
    size += sizeof(Val*) * keys_count;    // values pointers
    for (int i = 0; i < keys_count; i++)
      size += max_vals_per_key[i] * sizeof(Val);
    return size;
  }

 private:
  size_t keys_count_;
  Key* keys_;
  size_t* vals_per_key_;     // max size of circular array
  size_t* offsets_per_key_;  // one circular array per key
  Val** vals_;
};  // class PriorityQueue

};  // namespace linear
};  // namespace tools
};  // namespace neurox
