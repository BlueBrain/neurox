#pragma once

#include "neurox/neurox.h"

/**
 * @brief The Neuron class
 * Represents a neuron as its metadata and a tree-based orphology
 */
class Neuron
{
  public:
    hpx_t topBranch;		///> hpx address of the top compartment (soma)

    //neuron metadata
    int id;					///> neuron global id
};
