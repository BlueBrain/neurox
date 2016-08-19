#pragma once

#define DEBUG //Debug information at set-up time

//Core datatypes
#include "neurox/Neurox_hpx.h"
#include "neurox/input/InputParams.h"
#include "neurox/datatypes/Event.h"
#include "neurox/datatypes/VecPlayContinuous.h"
#include "neurox/datatypes/NetCon.h"
#include "neurox/datatypes/Mechanism.h"
#include "neurox/datatypes/Branch.h"
#include "neurox/datatypes/Neuron.h"

//Fixed-step Backward-Euler solver
#include "neurox/solver/BackwardEuler.h"
#include "neurox/solver/HinesSolver.h"

//CoreNeuron data loader
#include "neurox/input/coreneuron/Compartment.h"
#include "neurox/input/coreneuron/DataLoader.h"

//Miscellaneous
#include "neurox/misc/Statistics.h"

#define multiSpliX false
#define PARALLEL_MECHS_DEPENDENCY false

#define capacitance 3
#define IClamp 7
#define ProbAMPANMDA_EMS 137
#define ProbGABAAB_EMS 139

///NeuroX namespace: contains global information that is copied to all localities
namespace NeuroX
{
    extern int neuronsCount; 	///> total neurons count in the system
    extern hpx_t neuronsAddr; 	///> hpx address of the first position of the neurons array

    extern int mechanismsCount; ///> number of mechanisms
    extern NeuroX::Mechanism ** mechanisms; ///> array to all existing mechanisms
    extern int * mechanismsMap; ///>map of mechanisms offset in 'mechanisms' by 'mechanism type'

    extern Input::InputParams * inputParams; ///> Parameters parsed from command line

    extern hpx_action_t main;           ///> execution starting point (called via hpx_run)
    extern hpx_action_t clear;          ///> clears all memory utilised including neurons, branches and mechanisms information
    extern hpx_action_t setNeurons;     ///> Initialized neurons and neuronsAddr global vars
    extern hpx_action_t setInputParams; ///> Initializes InputParams
    extern hpx_action_t setMechanisms;	///> Initializes Mechanisms

    hpx_t getNeuronAddr(int i); ///> get HPX address of i-th neuron
    NeuroX::Mechanism * getMechanismFromType(int type); ///> returns mechanisms of type 'type'
    void registerHpxActions(); ///> Register all HPX actions

    static int main_handler(char **argv, size_t argc);
    static int clear_handler();
    static int setNeurons_handler(const int nargs, const void *args[], const size_t []);
    static int setInputParams_handler(const Input::InputParams * data, const size_t);
    static int setMechanisms_handler (const int nargs, const void *args[], const size_t[]);
} ;
