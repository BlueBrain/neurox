#pragma once

#include <utility>  //std::pair
#include "neurox/neurox.h"

namespace neurox {

/// Enumerator of all (used) event types
enum EventType
{
  kDiscreteEvent=0,
  //kTstopEvent=1,
  kNetCon=2,
  //kSelfEvent=3,
  //kPreSyn=4,
  //kNetParEvent=7,
  //kInputPreSyn=20,
  kVecPlayContinuous=-4
};

class Branch;

/**
 * @brief The Event virtual class
 * Simillar to DiscreteEvent in Coreneuron (netcon.h)
 * Stores the function calls that all events must have;
 */
class Event {
 public:
  virtual void Deliver(floble_t t, Branch* branch) = 0;
  virtual EventType Type() { return EventType::kDiscreteEvent; }
};
}

typedef std::pair<floble_t, neurox::Event*> TimedEvent;
