#pragma once

#include <map>
#include <vector>
#include <cassert>
#include <cstring> //memcpy

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

  typedef struct KeyInfoStruct {
    Key key_; /*keep as 1st field, to make std::bsearch below work*/
    size_t size_;  // max size of circular array
    Val* vals_; //Values
  } KeyInfo;

  /// Constructor for key to single-values
  Map(std::map<Key, Val>& map1, unsigned char* buffer) {
    // needs to be called with placement-new where buffer*
    assert((void*)buffer == this);

    vals_count_=0;
    keys_count_ = map1.size();
    assert(keys_count_>0);
    size_t offset = sizeof(Map<Key, Val>) + sizeof(KeyInfo)* keys_count_;

    int k = 0;
    KeyInfo *ki=nullptr;
    for (auto key_val : map1) {
      ki = &keys_info_[k];
      ki->key_ = key_val.first;
      ki->size_ = 1;
      ki->vals_ = (Val*)&(buffer[offset]);
      memcpy(ki->vals_, &(key_val.second), sizeof(Val));
      offset += sizeof(Val);
      vals_count_++;
      k++;
    }

    // make sure keys are sorted and not identical
    for (size_t k=1; k<keys_count_; k++) {
      assert(keys_info_[k].key_ > keys_info_[k-1].key_);
      if (keys_info_[k].key_ <= keys_info_[k-1].key_)
        throw std::runtime_error(
            std::string("Keys for linear map are not sorted."));
    }
  }

  /// Constructor for key array of values
  Map(std::map<Key, std::vector<Val*> >& map1, unsigned char* buffer) {
    // needs to be called with placement-new where buffer*
    assert((void*)buffer == this);

    vals_count_=0;
    keys_count_ = map1.size();
    assert(keys_count_>0);
    size_t offset = sizeof(Map<Key, Val>) + sizeof(KeyInfo) * keys_count_;

    int k=0, j=0;
    KeyInfo *ki=nullptr;
    for (auto key_val : map1) {
      ki = &keys_info_[k];
      ki->key_ = key_val.first;
      ki->size_ = key_val.second.size();
      ki->vals_ = (Val*)&(buffer[offset]);
      j = 0;
      for (auto val_ptr : key_val.second)
        memcpy(&(ki->vals_[j++]), val_ptr, sizeof(Val));
      offset += sizeof(Val) * ki->size_;
      vals_count_+=ki->size_;
      k++;
    }

    // make sure keys are sorted and not identical
    for (size_t k=1; k<keys_count_; k++) {
      assert(keys_info_[k].key_ > keys_info_[k-1].key_);
      if (keys_info_[k].key_ <= keys_info_[k-1].key_)
        throw std::runtime_error(
            std::string("Keys for linear map are not sorted."));
    }
  }

  ~Map() {
    for (int i = 0; i < keys_count_; i++) delete keys_info_[i]->vals_;
    delete[] keys_info_;
  }

  /// Buffer size for single value per key
  static size_t Size(size_t keys_count) {
    size_t size = sizeof(Map<Key, Val>) + sizeof(KeyInfo)*keys_count;
    size += sizeof(Val) * keys_count;     // single-value per key
    return size;
  }

  /// Buffer size of multiple values per key
  static size_t Size(size_t keys_count, size_t* vals_per_key) {
    size_t size = sizeof(Map<Key, Val>) + sizeof(KeyInfo)*keys_count;
    for (int i = 0; i < keys_count; i++) 
      size += vals_per_key[i] * sizeof(Val);
    return size;
  }

  /// single-value At
  inline Val* At(Key key) {
    KeyInfo* ki = (KeyInfo*)std::bsearch((void*)&key, (void*)keys_info_, keys_count_,
                                      sizeof(KeyInfo), Map::CompareKeyInfoPtrs);
    if (ki == nullptr)
      throw std::runtime_error(
          std::string("Key not found in linear map: " + key));

    return ki->vals_; //i.e. &ki->vals_[0]
  }

  /// multiple-value At
  inline void At(Key key, size_t& count, Val*& vals) {
    KeyInfo* ki = (KeyInfo*)std::bsearch((void*)&key, (void*)keys_info_, keys_count_,
                                      sizeof(KeyInfo), Map::CompareKeyInfoPtrs);
    if (ki == nullptr)
      throw std::runtime_error(
          std::string("Key not found in linear map: " + key));

    count = ki->size_;
    vals = ki->vals_;
  }

  KeyInfo* KeysInfo() { return keys_info_; }
  size_t Count() { return keys_count_; }
  Key KeyAt(size_t i) { assert(i<keys_count_); return keys_info_[i].key_; }
  Val* Values() { return keys_info_[0].vals_;} //values are serialized
  size_t ValuesCount() {return vals_count_;}

 private:
  size_t keys_count_;
  size_t vals_count_;
  KeyInfo *keys_info_;

  static int CompareKeyInfoPtrs(const void* pa, const void* pb) {
    const KeyInfo *a = (const KeyInfo *)pa;
    const KeyInfo *b = (const KeyInfo *)pb;
    return a->key_ < b->key_ ? -1 : (a->key_ == b->key_ ? 0 : 1);
  }
};  // class Map

};  // namespace linear
};  // namespace tools
};  // namespace neurox
