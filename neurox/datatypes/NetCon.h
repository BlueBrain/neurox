#pragma once

#include "neurox/Neurox.h"

namespace Neurox
{

/**
 * @brief The NetCon class
 * Equivalent to Coreneuron's NetCon in netcon.h
 * Includes the synaptic information at the post-synaptic neuron
 */
class NetConX
{
  public:
    NetConX();
    NetConX(int mechType, int mechInstance, double delay, double * args, short int argsCount, char active);
    ~NetConX();

    int mechType;   ///> mechanism type associated with this synapse
    short int argsCount;  ///> size of variable args
    char active;          ///> decides whether NetCon is active (or not)
    double delay;            ///> synaptic (soma-bouton distance + transmitters release delay)
    double * args;        ///> net_receive arguments (equivalent to weights in Coreneuron's NetCon)
    int mechInstance;     ///> mechanism instance, from the mechanism type

  private:
};

}
