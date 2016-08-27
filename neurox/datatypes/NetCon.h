#pragma once

#include "neurox/neurox.h"

namespace neurox
{

/**
 * @brief The NetCon class
 * Equivalent to Coreneuron's NetCon in netcon.h
 * Includes the synaptic information at the post-synaptic neuron
 */
class NetConX : Event
{
  public:
    NetConX();
    NetConX(int mechType, int mechInstance, double delay, double * args, short int argsCount, bool active);
    ~NetConX();

    void deliver(double t, Branch* branch) override; //event method (inherited)

    int mechType;   ///> mechanism type associated with this synapse
    short int argsCount;  ///> size of variable args
    bool active;          ///> decides whether NetCon is active (or not)
    double delay;         ///> synaptic (soma-bouton distance + transmitters release delay)
    double * args;        ///> net_receive arguments (equivalent to weights in Coreneuron's NetCon)
    int mechInstance;     ///> mechanism instance, from the mechanism type

  private:
};

///temp wrapper for point process
struct PointProcInfo
{
    int nodeId;
    int mechType;
    int mechInstance;
    int instanceDataOffset;
    int size;
};

}
