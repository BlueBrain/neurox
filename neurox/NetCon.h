#pragma once

#include "neurox/neurox.h"

namespace neurox {

/**
 * @brief The NetCon class
 * Equivalent to Coreneuron's NetCon in netcon.h
 * Includes the synaptic information at the post-synaptic neuron
 */
class NetConX : Event {
 public:
  NetConX();
  NetConX(int mechType, offset_t mechInstance, floble_t delay,
          offset_t weightIndex, unsigned short weightsCount, bool active);
  ~NetConX();

  void Deliver(floble_t t, Branch* branch) override;  // event method
                                                      // (inherited)

  EventType Type() { return EventType::kNetCon; }

  int mechType;                 ///> mechanism type associated with this synapse
  unsigned short weightsCount;  ///> size of variable args
  bool active;                  ///> decides whether NetCon is active (or not)

  /// synaptic delay (soma-bouton distance + transmitters release delay)
  floble_t delay;

  ///  net_receive arguments (equivalent to weights in Coreneuron's NetCon)
  offset_t weightIndex;

  /// mechanism instance, from the mechanism type
  offset_t mechInstance;

 private:
};

/// temp wrapper for point process
struct PointProcInfo {
  offset_t nodeId;  // compartment id
  int mechType;
  offset_t mechInstance;
  offset_t instanceDataOffset;
  size_t size;
};
}
