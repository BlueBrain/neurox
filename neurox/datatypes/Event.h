#pragma once

#include "neurox/Neurox.h"

#define DiscreteEventType 0
#define TstopEventType 1
#define NetConType 2
#define SelfEventType 3
#define PreSynType 4
#define NetParEventType 7
#define InputPreSynType 20
#define VecPlayContinuousType -4 //TODO is 4 but is already taken?

namespace NeuroX
{

class Branch;

/**
 * @brief The Event virtual class
 * Simullar DiscreteEvent in Coreneuron (netcon.h)
 * Stores the function calls that all events must have;
 */
class Event {
public:
    virtual void deliver(double t, Branch* branch)=0;
    virtual int type() { return DiscreteEventType; }
};

}
