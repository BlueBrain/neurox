/*
# =============================================================================
# Copyright (c) 2015 - 2021-2022 Blue Brain Project/EPFL
#
# See top-level LICENSE file for details.
# =============================================================================.
*/
#pragma once

#include <algorithm>  // std::sort
#include <cstring>    //memcpy

namespace neurox {
namespace tools {
namespace linear {

/* Use-case specific: VecPlayContinuousX has no neuron
 * id, so we give it a very large one. vecplay_max_vals
 * should be 1 (?), but we set to 10 to keep initial
 * events that are to be delivered at negative times */
static const bool add_vecplay_continuousx_entry = false;
static const int vecplay_max_vals = 10;
static const int vecplay_event_id = 999999999;

/**
 * @brief The Linear (max fixed size) PriorityQueue class
 */
template <class Key, class Time, class Val>
class PriorityQueue {
 public:
  // needs to be called with placement-new where buffer*
  // is the start of the data structure
  PriorityQueue() = delete;

  typedef struct KeyInfoStruct {
    Key key_;             /*keep as 1st field, to make std::bsearch below work*/
    size_t max_size_;     // max size of circular array
    size_t offset_pop_;   // one circular array per key
    size_t offset_push_;  // one circular array per key
    Val* vals_;           // Values
  } KeyInfo;

  PriorityQueue(size_t keys_count, Key* keys, size_t* max_vals_per_key,
                unsigned char* buffer) {
    assert((void*)buffer == this);

    this->keys_count_ = keys_count + (add_vecplay_continuousx_entry ? 1 : 0);
    KeyInfo* ki = nullptr;

    size_t offset = sizeof(PriorityQueue<Key, Time, Val>);
    this->keys_info_ = (KeyInfo*)&(buffer[offset]);
    for (int k = 0; k < keys_count_; k++) {
      ki = &keys_info_[k];
      if (k == keys_count_ - 1 && add_vecplay_continuousx_entry) {
        ki->key_ = (Key)vecplay_event_id;
        ki->max_size_ = vecplay_max_vals;
      } else {
        ki->key_ = keys[k];
        ki->max_size_ = max_vals_per_key[k];
      }
      ki->offset_push_ = 0;
      ki->offset_pop_ = 0;
    }
    offset += sizeof(KeyInfo) * keys_count_;

    Val dummy_val;
    for (int k = 0; k < keys_count_; k++) {
      ki = &keys_info_[k];
      ki->vals_ = (Val*)&(buffer[offset]);
      for (int j = 0; j < ki->max_size_; j++) {
        std::memcpy(&(ki->vals_[j]), &dummy_val, sizeof(Val));
        offset += sizeof(Val);
      }
    }

#ifndef NDEBUG
    // make sure keys are sorted and not identical
    for (int k = 1; k < keys_count_; k++) {
      assert(keys_info_[k].key_ > keys_info_[k - 1].key_);
      if (keys_info_[k].key_ <= keys_info_[k - 1].key_)
        throw std::runtime_error(
            std::string("Keys for linear priority queue are not sorted."));
    }
#endif
  }

  ~PriorityQueue() {
    for (int k = 0; k < keys_count_; k++) delete[] keys_info_[k]->vals_;
    delete[] keys_info_;
  }

  static size_t Size(size_t keys_count, size_t* max_vals_per_key) {
    size_t size =
        sizeof(PriorityQueue<Key, Time, Val>) + sizeof(KeyInfo) * keys_count;
    for (int i = 0; i < keys_count; i++)
      size += max_vals_per_key[i] * sizeof(Val);

    if (add_vecplay_continuousx_entry)
      size += sizeof(KeyInfo) + vecplay_max_vals * sizeof(Val);
    return size;
  }

  void Push(Key key, Val timed_event) {
    KeyInfo* ki = (KeyInfo*)std::bsearch(
        (void*)&key, (void*)keys_info_, keys_count_, sizeof(KeyInfo),
        PriorityQueue<Key, Time, Val>::CompareKeyInfoPtrs);
#ifndef NDEBUG
    if (ki == nullptr)
      throw std::runtime_error(
          std::string("Key not found in linear priority queue: " + key));
    assert(ki->key_ == key);
#endif

    size_t& offset_push = ki->offset_push_;
    std::memcpy(&(ki->vals_[offset_push]), &timed_event, sizeof(Val));
    if (++offset_push == ki->max_size_) offset_push = 0;
  }

  inline void PopAllBeforeTime(Time t, std::vector<Val>& events) {
    events.clear();
    KeyInfo* ki = nullptr;
    for (int k = 0; k < keys_count_; k++) {
      ki = &keys_info_[k];
      size_t& offset_pop = ki->offset_pop_;
      while (offset_pop != ki->offset_push_ &&
             ki->vals_[offset_pop].first <= t) {
        events.push_back(ki->vals_[offset_pop]);
        if (++offset_pop == ki->max_size_) offset_pop = 0;
      }
    }
    std::sort(events.begin(), events.end());
  }

  inline bool Empty() {
    KeyInfo* ki = nullptr;
    for (int i = 0; i < keys_count_; i++) {
      ki = &keys_info_[i];
      if (ki->offset_push_ != ki->offset_pop_[i]) return false;
    }
    return true;
  }

  KeyInfo* KeysInfo() { return keys_info_; }
  size_t Count() { return keys_count_; }

 private:
  size_t keys_count_;
  KeyInfo* keys_info_;

  static int CompareKeyInfoPtrs(const void* pa, const void* pb) {
    const KeyInfo* a = (const KeyInfo*)pa;
    const KeyInfo* b = (const KeyInfo*)pb;
    return a->key_ < b->key_ ? -1 : (a->key_ == b->key_ ? 0 : 1);
  }
};  // class PriorityQueue

};  // namespace linear
};  // namespace tools
};  // namespace neurox
