#pragma once

#define DEBUG

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

///Neurox namespace: contains global information that is copied to all localities
namespace Neurox
{
    extern int neuronsCount; 	///> total neurons count in the system
    extern hpx_t neuronsAddr; 	///> hpx address of the first position of the neurons array

    extern int mechanismsCount; ///> number of mechanisms
    extern Neurox::Mechanism * mechanisms; ///> Unique information per mechanism type
    extern Input::InputParams * inputParams; ///> Parameters parsed from command line

    extern hpx_action_t main;               ///> execution starting point (called via hpx_run)
    extern hpx_action_t setNeurons;     ///> Initialized neurons and neuronsAddr global vars
    extern hpx_action_t setInputParams;	    ///> Initializes InputParams
    extern hpx_action_t setMechanisms;	    ///> Initializes Mechanisms

    inline static hpx_t getNeuronAddr(int i) {
        return hpx_addr_add(neuronsAddr, sizeof(Neuron)*i, sizeof(Neuron));
    } ///> returns hpx address for i-th neuron

    void registerHpxActions(); ///> Register all HPX actions

    static int main_handler( char **argv, size_t argc);
    static int setNeurons_handler(const int nargs, const void *args[], const size_t []); ///>handler of setNeuronsAddr
    static int setInputParams_handler(const Input::InputParams * data, const size_t); ///> handler for setInputParams
    static int setMechanisms_handler (const int nargs, const void *args[], const size_t[]); ///> handler for setMechanisms
} ;
