#pragma once

#include "neurox/neurox.h"

namespace neurox {
namespace tools {
namespace linear {

/**
 * @brief The Linear Map class
 */
template <class Key, class Val>
class Map {
 public:
  // needs to be called with placement-new where buffer*
  // is the start of the data structure
  Map(size_t keys_count, Key* keys, size_t* vals_per_key, Val* vals,
      char* buffer)
      : keys_count_(keys_count) {
    assert(buffer == this);

    // keys_
    size_t offset = sizeof(size_t);
    keys_ = (Key*)&(buffer[offset]);
    memcpy(keys_, keys, sizeof(Key) * keys_count);

    // count of values per key
    offset += sizeof(Key) * keys_count;
    vals_per_key_ = (size_t*)&(buffer[offset]);
    memcpy(vals_per_key_, vals_per_key, sizeof(size_t) * keys_count);

    // values array of pointers
    offset += sizeof(size_t) * keys_count;
    vals_ = (Val**)&(buffer[offset]);

    // values
    offset += sizeof(Val*) * keys_count;
    size_t vals_offset = 0;
    for (int i = 0; i < keys_count; i++) {
      vals_[i] = &buffer[offset];
      memcpy(vals_[i], &vals[vals_offset], sizeof(Val) * vals_per_key[i]);
      offset += vals_per_key[i] * sizeof(Val);
    }
  }

  size_t Size(size_t keys_count, size_t* vals_per_key) {
    size_t size = sizeof(Map<Key, Val>);
    size += sizeof(size_t) * keys_count;  // keys
    size += sizeof(size_t) * keys_count;  // vals per key
    size += sizeof(Val*) * keys_count;    // values pointers
    for (int i = 0; i < keys_count; i++) size += vals_per_key[i] * sizeof(Val);
    return size;
  }

  inline size_t GetIndex(Key key) { return -1; }

  inline void At(Key key, size_t& count, Val*& vals) {
    size_t i = GetIndex(key);
    count = vals_per_key_[i];
    vals = vals_[i];
  }

 private:
  size_t keys_count_;
  Key* keys_;
  size_t* vals_per_key_;
  Val** vals_;
};  // class Map

};  // namespace linear
};  // namespace tools
};  // namespace neurox
