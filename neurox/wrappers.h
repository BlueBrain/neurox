#pragma once

namespace neurox {
namespace wrappers {

/// Memory pinning from an hpx memorry address to a local pointer
template <typename T>
inline hpx_t MemoryPin(T) {
  hpx_t target = hpx_thread_current_target();
  T* local = NULL;
  return hpx_gas_try_pin(target, (void**)&local);
}

/// Memory UNpinning from an hpx memorry address to a local pointer
inline hpx_status_t MemoryUnpin(hpx_t target) {
  hpx_gas_unpin(target);
  return HPX_SUCCESS;
}

/// Memory UNpinning from an hpx memorry address to a local pointer
/// with the return of a value
template <typename T>
inline hpx_status_t MemoryUnpin(hpx_t target, const T& var) {
  hpx_gas_unpin(target);
  return HPX_THREAD_CONTINUE(var);
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

/// call action (with arguments) on all localities
template <typename... Args>
inline hpx_status_t CallAllLocalities(hpx_action_t f, Args... args) {
  return hpx_bcast_rsync(f, args...);
}

/// calls method (with arguments) on all neurons in neurox::neurons
template <typename... Args>
inline hpx_status_t CallAllNeurons(hpx_action_t f, Args... args) {
  hpx_t lco = hpx_lco_and_new(neurox::neurons_count);
  int e = HPX_SUCCESS;
  int n = wrappers::CountArgs(args...);
  for (size_t i = 0; i < neurox::neurons_count; i++)
    e += _hpx_call(neurox::neurons[i], f, lco, n, args...);
  hpx_lco_wait_reset(lco);
  hpx_lco_delete_sync(lco);
  return e;
}

/// register hpx-action and handlers for zero-variables action
static void RegisterZeroVarAction(hpx_action_t& action, int (*handler)(void)) {
  HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE, action, handler);
}

/// register hpx-action and handlers for action with single variable
template <typename T>
static void RegisterSingleVarAction(hpx_action_t& action,
                                    int (*handler)(const T*, const size_t),
                                    bool compressed = false) {
  HPX_REGISTER_ACTION(HPX_DEFAULT,
                      HPX_MARSHALLED | (compressed ? HPX_COMPRESSED : 0),
                      action, handler, HPX_POINTER, HPX_SIZE_T);
}

/// register hpx-action and handlers for action with multiples variables
static void RegisterMultipleVarAction(hpx_action_t& action,
                                      int (*handler)(const int, const void* [],
                                                     const size_t[]),
                                      bool compressed = false) {
  HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED | HPX_VECTORED |
                                       (compressed ? HPX_COMPRESSED : 0),
                      action, handler, HPX_INT, HPX_POINTER, HPX_POINTER);
}

/// register hpx-action and handlers for an AllReduce init action
template <typename T>
static void RegisterAllReduceInitAction(hpx_action_t& action,
                                        void (*handler)(T*, const size_t)) {
  HPX_REGISTER_ACTION(HPX_FUNCTION, 0, action, handler);
}

/// register hpx-action and handlers for an AllReduce init action
template <typename T>
static void RegisterAllReduceReduceAction(hpx_action_t& action,
                                          void (*handler)(T*, const T*,
                                                          const size_t)) {
  HPX_REGISTER_ACTION(HPX_FUNCTION, 0, action, handler);
}

/// get running thread Id
inline int MyThreadId() { return hpx_thread_get_tls_id(); }

/// get current locality Id on the network
inline int MyRank() { return hpx_get_my_rank(); }

/// get number of localities in the network
inline int NumRanks() { return hpx_get_num_ranks(); }

/// get number of threads in current locality
inline int NumThreads() { return hpx_get_num_threads(); }

};  // namespace wrappers;
};  // namespare neurox;
