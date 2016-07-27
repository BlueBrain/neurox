#pragma once

//Core datatypes
#include "neurox/Neurox_hpx.h"
#include "neurox/datatypes/Capacitance.h"
#include "neurox/input/InputParams.h"
#include "neurox/datatypes/NetCon.h"
#include "neurox/datatypes/Spike.h"
#include "neurox/datatypes/Mechanism.h"
#include "neurox/datatypes/Branch.h"
#include "neurox/datatypes/Neuron.h"

//Solvers
#include "neurox/solver/BackwardEuler.h"

//CoreNeuron data loader
#include "neurox/input/coreneuron/Compartment.h"
#include "neurox/input/coreneuron/DataLoader.h"

namespace Neurox
{
    extern int neuronsCount; 	///> total neurons count in the system
    extern hpx_t neuronsAddr; 	///> hpx address of the first position of the neurons array

    extern int mechanismsCount; ///> number of mechanisms
    extern Neurox::Mechanism * mechanisms; ///> Unique information per mechanism type

    extern Input::InputParams * inputParams; ///> Parameters parsed from command line (TODO should go away at some point)

    inline static hpx_t getNeuronAddr(int i) {
        return hpx_addr_add(neuronsAddr, sizeof(Neuron)*i, sizeof(Neuron));
    } ///> returns hpx address for i-th neuron

    static hpx_action_t setInputParams;	///> Initializes InputParams
    static hpx_action_t setMechanisms;	///> Initializes Mechanisms

    int setInputParams_handler(const Input::InputParams * inputParams, const size_t size); ///handler for setInputParams
    int setMechanisms_handler(const Mechanism * mechanisms, const size_t size); ///handler for setMechanisms
} ;
