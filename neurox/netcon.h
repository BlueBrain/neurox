#pragma once

#include "neurox/neurox.h"

namespace neurox {

/**
 * @brief The NetCon class
 * Equivalent to Coreneuron's NetCon in netcon.h
 * Includes the synaptic information at the post-synaptic neuron
 */
class NetconX : Event {
 public:
  NetconX();
  NetconX(int mech_type, offset_t mech_instance, floble_t delay_,
          offset_t weight_index, unsigned short weights_count, bool active_);
  ~NetconX();

  void Deliver(floble_t t, Branch* branch) override;  // event method
                                                      // (inherited)

  EventTypes Type() { return EventTypes::kNetCon; }

  int mech_type_;  ///> mechanism type associated with this synapse
  unsigned short weights_count_;  ///> size of variable args
  bool active_;                   ///> decides whether NetCon is active (or not)

  /// synaptic delay (soma-bouton distance + transmitters release delay)
  floble_t delay_;

  ///  net_receive arguments (equivalent to weights in Coreneuron's NetCon)
  offset_t weight_index_;

  /// mechanism instance, from the mechanism type
  offset_t mech_instance_;

 private:
};

/// temp wrapper for point process
struct PointProcInfo {
  offset_t node_id;  // compartment id
  int mech_type;
  offset_t mech_instance;
  offset_t instance_data_offset;
  size_t size;
};
}
