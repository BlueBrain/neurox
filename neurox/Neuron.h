#pragma once

#include <deque>
#include "libhpx/libhpx.h"
#include "neurox/neurox.h"

namespace neurox {

/**
 * @brief The Neuron class
 * Represents the soma structure and outgoing spikes network of a neuron
 */
class Neuron {
 public:
  Neuron() = delete;
  ~Neuron();

  Neuron(neuron_id_t neuronId, floble_t APthreshold);

  neuron_id_t gid;     ///> neuron global id
  floble_t threshold;  ///> Action Potential threshold (PreSyn->_threshold)
  floble_t refractoryPeriod;  ///> refractory period

  void setupTreeMatrixMinimal();  ///>set_tree_matrix_minimal

  /// checks if AP threshold has been reached and if spikes can be transmitted
  /// (PreSynHelper->flag)
  bool CheckAPthresholdAndTransmissionFlag(floble_t v);

  /// the outgoing synapse:
  class Synapse {
   public:
    Synapse() = delete;
    Synapse(hpx_t branchAddr, floble_t minDelay, hpx_t topBranchAddr = HPX_NULL,
            int destinationGid = -1);
    ~Synapse();
    hpx_t branchAddr;     ///> address of destination
    hpx_t topBranchAddr;  ///> addres of top-branch (soma) of destination neuron
#ifndef NDEBUG
    int destinationGid;
#endif
    ///  next time this post-syn neuron needs to be informed of my actual time
    floble_t nextNotificationTime;
    /// interval  of notification in case of no spykes (fastest Netcon from
    /// current neuron to dependant-neuron)
    floble_t minDelay;
    hpx_t previousSpikeLco;  ///>lco controlling spikes delivery
  };

  hpx_t SendSpikes(floble_t t);  ///> fires AP, returns LCO for sent synapses
  void AddSynapse(
      Synapse* target);       ///> add hpx address of post-synaptic branch
  size_t GetSynapsesCount();  ///> get size of vector synapse

  algorithms::AlgorithmMetaData* algorithmMetaData;

  std::vector<Synapse*> synapses;  ///> out-going synapse information

 private:
  // the outgoing neuron connections:
  hpx_t synapsesMutex;  ///> protects synapses

  ///  PreSynHelper* psh->flag (whether spikes for a given AP have been sent)
  bool synapsesTransmissionFlag;
};  // Neuron
};  // namespace neurox
