#pragma once

namespace neurox {
namespace wrappers {

/// Memory pinning from an hpx memorry address to a local pointer
template <typename T>
bool MemoryPin(T *&local, hpx_t &target) {
  target = hpx_thread_current_target();
  return hpx_gas_try_pin(target, (void **)&local);
}

/// Memory UNpinning from an hpx memorry address to a local pointer
template <typename T>
hpx_status_t MemoryUnpin(hpx_t target) {
  hpx_gas_unpin(target);
  return HPX_SUCCESS;
}

/// call action (with arguments) on all localities
template <typename... Args>
hpx_status_t CallAllLocalities(hpx_action_t f, Args... args) {
  return hpx_bcast_rsync(f, args...);
}

/// count the number of arguments
template <typename... ArgTypes>
inline int CountArgs(ArgTypes... args);

template <typename T, typename... ArgTypes>
inline int CountArgs(T t, ArgTypes... args) {
  return 1 + CountArgs(args...);
}
template <>
inline int CountArgs() {
  return 0;
}

/// calls method (with arguments) on all neurons in neurox::neurons
template <typename... Args>
hpx_status_t CallAllNeurons(hpx_action_t f, Args... args) {
  hpx_t lco = hpx_lco_and_new(neurox::neurons_count);
  int e = HPX_SUCCESS;
  int n = neurox::wrappers::CountArgs(args...);
  for (size_t i = 0; i < neurox::neurons_count; i++)
    e += _hpx_call(neurox::neurons[i], f, lco, n, args...);
  hpx_lco_wait_reset(lco);
  hpx_lco_delete_sync(lco);
  return e;
}

}; //namespace wrappers;
}; //namespare neurox;
