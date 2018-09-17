#pragma once

#include <map>
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
  Map() = delete;

  Map(std::map<Key, std::vector<Val*> >& map1, unsigned char* buffer) {
    // needs to be called with placement-new where buffer*
    // is the start of the data structure
    assert(buffer == this);
    size_t offset = 0;

    keys_count_ = map1.size();
    offset += sizeof(size_t);
    keys_ = (Key*)&(buffer[offset]);
    offset += sizeof(Key) * keys_count_;
    vals_per_key_ = (size_t*)&(buffer[offset]);
    offset += sizeof(size_t) * keys_count_;
    vals_ = (Val**)&(buffer[offset]);
    offset += sizeof(Val*) * keys_count_;

    int i = 0;
    for (auto map_it : map1) {
      keys_[i] = map_it.first;
      vals_per_key_[i] = map_it.second.size();
      vals_[i] = (Val*)&(buffer[offset]);
      int j = 0;
      for (auto vec_it : map_it.second)
        memcpy(&(vals_[i][j++]), vec_it, sizeof(Val));
      i++;
    }
  }

  static size_t Size(size_t keys_count, size_t* vals_per_key) {
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
