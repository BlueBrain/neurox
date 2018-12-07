#pragma once

#include <cstring>	    //memcpy
#include <algorithm>    // std::sort

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

  PriorityQueue(size_t keys_count, Key* keys, size_t* max_vals_per_key,
                unsigned char* buffer) {
    assert((void*)buffer == this);

    size_t * delete_ptr = nullptr;
    if (add_vecplay_continuousx_entry) {
      Key* keys2 = new Key[keys_count + 1];
      size_t* max_vals_per_key2 = new size_t[keys_count + 1];
      delete_ptr = max_vals_per_key2; //fix mem leak
      std::memcpy(keys2, keys, sizeof(Key) * keys_count);
      std::memcpy(max_vals_per_key2, max_vals_per_key, sizeof(size_t) * keys_count);
      keys2[keys_count] = (Key)vecplay_event_id;
      max_vals_per_key2[keys_count] = vecplay_max_vals;
      keys = keys2;
      max_vals_per_key = max_vals_per_key2;
      keys_count++;
    }

    keys_count_ = keys_count;
    size_t offset = sizeof(PriorityQueue<Key, Time, Val>);

    // keys_
    keys_ = (Key*)&(buffer[offset]);
    std::memcpy(keys_, keys, sizeof(Key) * keys_count);
    offset += sizeof(Key) * keys_count;

    // count of values per key (size of circular array)
    max_vals_per_key_ = (size_t*)&(buffer[offset]);
    std::memcpy(max_vals_per_key_, max_vals_per_key, sizeof(size_t) * keys_count);
    offset += sizeof(size_t) * keys_count;

    // current offset of each circular array
    offsets_pop_ = (size_t*)&(buffer[offset]);
    for (int i = 0; i < keys_count; i++) offsets_pop_[i] = 0;
    offset += sizeof(size_t) * keys_count;

    offsets_push_ = (size_t*)&(buffer[offset]);
    for (int i = 0; i < keys_count; i++) offsets_push_[i] = 0;
    offset += sizeof(size_t) * keys_count;

    // values array of pointers
    vals_ = (Val**)&(buffer[offset]);
    offset += sizeof(Val*) * keys_count;

    // values
    for (int i = 0; i < keys_count; i++) {
      vals_[i] = (Val*)&(buffer[offset]);
      Val dummy_val;
      for (int j = 0; j < max_vals_per_key[i]; j++) {
        std::memcpy(&(vals_[i][j]), &dummy_val, sizeof(Val));
        offset += sizeof(Val);
      }
    }

    // make sure keys are sorted and not identical
    for (int i = 1; i < keys_count_; i++) {
      assert(keys_[i] > keys_[i - 1]);
      if (keys_[i] <= keys_[i - 1])
        throw std::runtime_error(
            std::string("Keys for linear priority queue are not sorted."));
    }

    if (delete_ptr)
      delete [] delete_ptr;
  }

  ~PriorityQueue() {
    delete[] keys_;
    delete[] max_vals_per_key_;
    delete[] offsets_push_;
    delete[] offsets_pop_;
    for (int i = 0; i < keys_count_; i++) delete[] vals_[i];
    delete[] vals_;
  }

  static size_t Size(size_t keys_count, size_t* max_vals_per_key) {
    size_t * delete_ptr = nullptr;
    if (add_vecplay_continuousx_entry) {
      size_t* max_vals_per_key2 = new size_t[keys_count + 1];
      delete_ptr = max_vals_per_key2; //mem leak
      std::memcpy(max_vals_per_key2, max_vals_per_key, sizeof(size_t) * keys_count);
      max_vals_per_key2[keys_count] = vecplay_max_vals;
      max_vals_per_key = max_vals_per_key2;
      keys_count++;
    }

    size_t size = sizeof(PriorityQueue<Key, Time, Val>);
    size += sizeof(Key) * keys_count;     // keys
    size += sizeof(size_t) * keys_count;  // vals per key
    size += sizeof(size_t) * keys_count;  // pop offsets per key
    size += sizeof(size_t) * keys_count;  // push offsets per key
    size += sizeof(Val*) * keys_count;    // values pointers
    for (int i = 0; i < keys_count; i++)
      size += max_vals_per_key[i] * sizeof(Val);

    if (delete_ptr)
      delete [] delete_ptr;
    return size;
  }

  void Push(Key key, Val timed_event) {
    Key* key_ptr =
        (Key*)std::bsearch((void*)&key, (void*)keys_, keys_count_, sizeof(Key),
                           PriorityQueue<Key, Time, Val>::CompareKeyPtrs);
    if (key_ptr == nullptr)
      throw std::runtime_error(
          std::string("Key not found in linear priority queue: " + key));
    assert(*key_ptr == key);

    const size_t k = key_ptr - keys_;
    const size_t max_vals = max_vals_per_key_[k];
    size_t& offset_push = offsets_push_[k];
    std::memcpy(&(vals_[k][offset_push]), &timed_event, sizeof(Val));
    if (++offset_push == max_vals) offset_push = 0;
  }

  void PopAllBeforeTime(Time t, std::vector<Val>& events) {
    events.clear();
    for (int k = 0; k < keys_count_; k++) {
      const size_t max_vals = max_vals_per_key_[k];
      const size_t offset_push = offsets_push_[k];
      size_t& offset_pop = offsets_pop_[k];
      while (offset_pop != offset_push && vals_[k][offset_pop].first <= t) {
        events.push_back(vals_[k][offset_pop]);
        if (++offset_pop == max_vals) offset_pop = 0;
      }
    }
    std::sort(events.begin(), events.end());
  }

  bool Empty() {
    for (int i = 0; i < keys_count_; i++)
      if (offsets_push_[i] > offsets_pop_[i]) return false;
    return true;
  }

  Key* Keys() { return keys_; }
  size_t Count() { return keys_count_; }

 private:
  static int CompareKeyPtrs(const void* pa, const void* pb) {
    const Key a = *(const Key*)pa;
    const Key b = *(const Key*)pb;
    return a < b ? -1 : (a == b ? 0 : 1);
  }
  size_t keys_count_;
  Key* keys_;
  size_t* max_vals_per_key_;  // max size of circular array
  size_t* offsets_pop_;       // one circular array per key
  size_t* offsets_push_;      // one circular array per key
  Val** vals_;
};  // class PriorityQueue

};  // namespace linear
};  // namespace tools
};  // namespace neurox
