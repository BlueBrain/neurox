#pragma once

#include "neurox/neurox.h"

#include <deque>

namespace neurox {

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
    Synapse(hpx_t branch_addr_, floble_t min_delay_,
            hpx_t top_branch_addr_ = HPX_NULL, int destination_gid_ = -1);
    ~Synapse();
    hpx_t branch_addr_;      ///> address of destination
    hpx_t top_branch_addr_;  ///> addres of top-branch (soma) of destination
                             /// neuron
#ifndef NDEBUG
    int destination_gid_;
#endif
    ///  next time this post-syn neuron needs to be informed of my actual time
    floble_t next_notification_time_;
    /// interval  of notification in case of no spykes (fastest Netcon from
    /// current neuron to dependant-neuron)
    floble_t min_delay_;
    hpx_t previous_spike_lco_;  ///>lco controlling spikes delivery
  };

  hpx_t SendSpikes(floble_t t);  ///> fires AP, returns LCO for sent synapses
  void AddSynapse(
      Synapse* target);       ///> add hpx address of post-synaptic branch
  size_t GetSynapsesCount();  ///> get size of vector synapse

  algorithms::AlgorithmMetadata* algorithm_metadata_;

  std::vector<Synapse*> synapses_;  ///> out-going synapse information

 private:
  // the outgoing neuron connections:
  hpx_t synapses_mutex_;  ///> protects synapses

  ///  PreSynHelper* psh->flag (whether spikes for a given AP have been sent)
  bool synapses_transmission_flag_;
};  // Neuron
};  // namespace neurox
