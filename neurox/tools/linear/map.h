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
    assert((void*)buffer == this);
    size_t offset = 0;

    keys_count_ = map1.size();
    offset += sizeof(Map<Key, Val>);
    keys_ = (Key*)&(buffer[offset]);
    offset += sizeof(Key) * keys_count_;
    vals_per_key_ = (size_t*)&(buffer[offset]);
    offset += sizeof(size_t) * keys_count_;
    vals_ = (Val**)&(buffer[offset]);
    offset += sizeof(Val*) * keys_count_;

    int i = 0;
    for (auto key_vals : map1) {
      keys_[i] = key_vals.first;
      vals_per_key_[i] = key_vals.second.size();
      vals_[i] = (Val*)&(buffer[offset]);
      int j = 0;
      for (auto val_ptr : key_vals.second)
        memcpy(&(vals_[i][j++]), val_ptr, sizeof(Val));
      offset += vals_per_key_[i] * sizeof(Val);
      i++;
    }
    assert(offset == Size(keys_count_, vals_per_key_));

    // make sure keys are sorted and not identical
    for (int i = 1; i < keys_count_; i++) {
      assert(keys_[i] > keys_[i - 1]);
      if (keys_[i] <= keys_[i - 1])
        throw std::runtime_error(
            std::string("Keys for linear map are not sorted."));
    }
  }

  ~Map() {
    delete[] keys_;
    delete[] vals_per_key_;
    for (int i = 0; i < keys_count_; i++) delete[] vals_[i];
    delete[] vals_;
  }

  static size_t Size(size_t keys_count, size_t* vals_per_key) {
    size_t size = sizeof(Map<Key, Val>);
    size += sizeof(Key) * keys_count;     // keys
    size += sizeof(size_t) * keys_count;  // vals per key
    size += sizeof(Val*) * keys_count;    // values pointers
    for (int i = 0; i < keys_count; i++) size += vals_per_key[i] * sizeof(Val);
    return size;
  }

  inline size_t GetIndex(Key key) { return -1; }

  inline void At(Key key, size_t& count, Val*& vals) {
    Key* k = (Key*)std::bsearch((void*)&key, (void*)keys_, keys_count_,
                                sizeof(Key), Map::CompareKeyPtrs);
    if (k == nullptr)
      throw std::runtime_error(
          std::string("Key not found in linear map: " + key));
    assert(*k == key);
    const size_t pos = k - keys_;
    count = vals_per_key_[pos];
    vals = vals_[pos];
  }

  Key* Keys() { return keys_; }
  size_t Count() { return keys_count_; }

 private:
  static int CompareKeyPtrs(const void* pa, const void* pb) {
    const Key a = *(const Key*)pa;
    const Key b = *(const Key*)pb;
    return a<b ? -1 : (a==b ? 0 : 1);
  }
  size_t keys_count_;
  Key* keys_;
  size_t* vals_per_key_;
  Val** vals_;
};  // class Map

};  // namespace linear
};  // namespace tools
};  // namespace neurox
