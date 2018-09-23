#pragma once

#define ADD_VECPLAY_CONTINUOUSX_ENTRY

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
                unsigned char* buffer) {
    assert((void*)buffer == this);

#ifdef ADD_VECPLAY_CONTINUOUSX_ENTRY
    Key* keys2 = new Key[keys_count + 1];
    size_t* max_vals_per_key2 = new size_t[keys_count + 1];
    memcpy(&keys2[1], keys, sizeof(Key) * keys_count);
    memcpy(&max_vals_per_key2[1], max_vals_per_key,
           sizeof(size_t) * keys_count);
    keys2[0] = (Key)-999999999;
    max_vals_per_key2[0] = 1;
    delete[] keys;
    delete[] max_vals_per_key;
    keys = keys2;
    max_vals_per_key = max_vals_per_key2;
    keys_count++;
#endif

    keys_count_ = keys_count;

    size_t offset = sizeof(PriorityQueue<Key, Val>);

    // keys_
    keys_ = (Key*)&(buffer[offset]);
    memcpy(keys_, keys, sizeof(Key) * keys_count);
    offset += sizeof(Key) * keys_count;

    // count of values per key (size of circular array)
    vals_per_key_ = (size_t*)&(buffer[offset]);
    memcpy(vals_per_key_, max_vals_per_key, sizeof(size_t) * keys_count);
    offset += sizeof(size_t) * keys_count;

    // current offset of each circular array
    offsets_push_ = (size_t*)&(buffer[offset]);
    for (int i = 0; i < keys_count; i++) offsets_push_[i] = 0;
    offset += sizeof(size_t) * keys_count;

    offsets_pop_ = (size_t*)&(buffer[offset]);
    for (int i = 0; i < keys_count; i++) offsets_pop_[i] = 0;
    offset += sizeof(size_t) * keys_count;

    // values array of pointers
    vals_ = (Val**)&(buffer[offset]);
    offset += sizeof(Val*) * keys_count;

    // values
    for (int i = 0; i < keys_count; i++) {
      vals_[i] = (Val*)&(buffer[offset]);
      Val dummy_val;
      for (int j = 0; j < max_vals_per_key[i]; j++) {
        memcpy(&(vals_[i][j]), &dummy_val, sizeof(Val));
        offset += sizeof(Val);
      }
    }
  }

  ~PriorityQueue() {
    delete[] keys_;
    delete[] vals_per_key_;
    delete[] offsets_push_;
    delete[] offsets_pop_;
    for (int i = 0; i < keys_count_; i++) delete[] vals_[i];
    delete[] vals_;
  }

  static size_t Size(size_t keys_count, size_t* max_vals_per_key) {
#ifdef ADD_VECPLAY_CONTINUOUSX_ENTRY
    size_t* max_vals_per_key2 = new size_t[keys_count + 1];
    memcpy(&max_vals_per_key2[1], max_vals_per_key,
           sizeof(size_t) * keys_count);
    max_vals_per_key2[0] = 1;
    delete[] max_vals_per_key;
    max_vals_per_key = max_vals_per_key2;
    keys_count++;
#endif
    size_t size = sizeof(PriorityQueue<Key, Val>);
    size += sizeof(Key) * keys_count;     // keys
    size += sizeof(size_t) * keys_count;  // vals per key
    size += sizeof(size_t) * keys_count;  // insert offsets per key
    size += sizeof(size_t) * keys_count;  // get offsets per key
    size += sizeof(Val*) * keys_count;    // values pointers
    for (int i = 0; i < keys_count; i++)
      size += max_vals_per_key[i] * sizeof(Val);
    return size;
  }

  void Push(Key key, Val timed_event) {
    for (int i = 0; i < keys_count_; i++)
      if (keys_[i] == key) {
        size_t& offset = offsets_push_[i];
        size_t count = vals_per_key_[i];
        memcpy(&(vals_[i][offset]), &timed_event, sizeof(Val));
        offset++;
        if (offset == count) offset = 0;
        break;
      }
    assert(0);
  }

  void PopAllBeforeTime(floble_t t, std::vector<Val>& events) {
    for (int i = 0; i < keys_count_; i++) {
      size_t& offset_pop = offsets_pop_[i];
      size_t& offset_push = offsets_push_[i];
      while (offset_pop < offset_push && vals_[i][offset_pop].first <= t) {
        events.push_back(vals_[i][offset_pop]);
        offset_pop++;
      }
    }
    std::sort(events.begin(), events.end());
  }

  Val* Pop(Key key) {
    for (int i = 0; i < keys_count_; i++)
      if (keys_[i] == key) {
        size_t& offset = offsets_pop_[i];

        // no event
        if (offset == offsets_push_[i]) return nullptr;

        // return value
        Val* ret_val = &(vals_[i][offset]);

        // append read offset
        size_t count = vals_per_key_[i];
        offset++;
        if (offset == count) offset = 0;

        return ret_val;
      }
    assert(0);
  }

  bool Empty() {
    for (int i = 0; i < keys_count_; i++)
      if (offsets_push_[i] > offsets_pop_[i]) return true;
    return true;
  }

 private:
  size_t keys_count_;
  Key* keys_;
  size_t* vals_per_key_;  // max size of circular array
  size_t* offsets_push_;  // one circular array per key
  size_t* offsets_pop_;   // one circular array per key
  Val** vals_;
};  // class PriorityQueue

};  // namespace linear
};  // namespace tools
};  // namespace neurox
