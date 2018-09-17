#pragma once

#include "neurox/neurox.h"

#include <deque>
#include <set>

using namespace tools;

namespace neurox {

/// forward declarations
namespace synchronizers {
class SynchronizerNeuronInfo;
};

/**
 * @brief The Neuron class
 * Represents the soma structure and outgoing spikes network of a neuron
 */
class Neuron {
 public:
  Neuron() = delete;
  ~Neuron();

  Neuron(neuron_id_t neuron_id, floble_t ap_threshold);

  neuron_id_t gid_;     ///> neuron global id
  floble_t threshold_;  ///> Action Potential threshold (PreSyn->_threshold)
  floble_t refractory_period_;  ///> refractory period

  void SetupTreeMatrixMinimal();  ///>set_tree_matrix_minimal

  /// checks if AP threshold has been reached and if spikes can be transmitted
  /// (PreSynHelper->flag)
  bool CheckAPthresholdAndTransmissionFlag(floble_t v);

  /// the outgoing synapse:
  class Synapse {
   public:
    Synapse() = delete;
    Synapse(hpx_t branch_addr, floble_t min_delay,
            hpx_t top_branch_addr = HPX_NULL, int destination_gid = -1);
    ~Synapse();

    /// address of destination branch or locality
    hpx_t branch_addr_;

    /// address of top-branch (soma) or locality of destination neuron
    hpx_t soma_or_locality_addr_;

#ifndef NDEBUG
    int destination_gid_;
#endif
    ///  next time this post-syn neuron needs to be informed of my actual time
    floble_t next_notification_time_;

    /// interval  of notification in case of no spykes (fastest Netcon from
    /// current neuron to dependant-neuron or dependant-locality)
    floble_t min_delay_;
    hpx_t previous_spike_lco_;  ///>lco controlling spikes delivery
  };

  /// fires AP, returns LCO for sent synapses
  hpx_t SendSpikes(floble_t t);

  /// returns whether progress of neuron is handled by a scheduler
  inline bool HasScheduler() { return synchronizer_step_trigger_ != HPX_NULL; }

  /// add hpx address of post-synaptic branch
  void AddSynapse(Synapse*);

  /// copy data from synapses_ to synapses_linear_
  void LinearizeSynapses();

  /// get size of vector synapse
  size_t GetSynapsesCount();

  /// the outgoing neuron connections:
  std::vector<Synapse*> synapses_;

  /// linear data container of synapses
  linear::Vector<Synapse>* synapses_linear_;

  /// Synchronizer-dependent metadata
  synchronizers::SynchronizerNeuronInfo* synchronizer_neuron_info_;

  /// and-gate (trigger) for time-based synchronization of stepping
  hpx_t synchronizer_step_trigger_;

 private:
  hpx_t synapses_mutex_;  ///> mutex protecting synapses

  ///  PreSynHelper* psh->flag (whether spikes for a given AP have been sent)
  bool synapses_transmission_flag_;

  unsigned char* synapses_linear_buffer_;
};  // Neuron
};  // namespace neurox
