#pragma once

//Coreneuron core datatypes
//#define NET_RECEIVE_BUFFERING 0 //disables buffering on net_receive functions. TODO past on CMake, should be put here instead
#include "coreneuron/coreneuron.h"

//Core datatypes
#include "neurox/hpx.h"
#include "neurox/Event.h"
#include "neurox/VecPlayContinuous.h"
#include "neurox/NetCon.h"
#include "neurox/Mechanism.h"
#include "neurox/Branch.h"
#include "neurox/Neuron.h"

//Fixed-step Backward-Euler solver
#include "neurox/input/InputParams.h"
#include "neurox/solver/HinesSolver.h"

//CoreNeuron Input
#include "neurox/input/coreneuron/Compartment.h"
#include "neurox/input/coreneuron/DataLoader.h"
#include "neurox/input/coreneuron/Debugger.h"

//Miscellaneous
#include "neurox/misc/Statistics.h"

#define IClamp 7
#define ProbAMPANMDA_EMS 137
#define ProbGABAAB_EMS 139

//#define NDEBUG //(if undefined, compares data output with coreneuron and enables assertions)
//#define PRINT_EVENT
#define PRINT_TIME_DEPENDENCY

///neurox namespace: contains global information that is copied to all localities
namespace neurox
{
    extern int neuronsCount; 	///> total neurons count in the system
    extern hpx_t neuronsAddr; 	///> hpx address of the neurons array

    extern int mechanismsCount; ///> number of mechanisms
    extern neurox::Mechanism ** mechanisms; ///> array to all existing mechanisms
    extern int * mechanismsMap; ///>map of mechanisms offset in 'mechanisms' by 'mechanism type'

    extern Input::InputParams *  inputParams; ///> Parameters parsed from command line

    extern hpx_action_t main;            ///> execution starting point (called via hpx_run)
    extern hpx_action_t clear;           ///> clears all memory utilised including neurons, branches and mechanisms information
    extern hpx_action_t setNeurons;      ///> Initialized neurons and neuronsAddr global vars
    extern hpx_action_t setInputParams;  ///> Initializes InputParams
    extern hpx_action_t setMechanisms;	 ///> Initializes Mechanisms

    hpx_t getNeuronAddr(int i);      ///> get HPX address of i-th neuron
    Mechanism * getMechanismFromType(int type); ///> returns mechanisms of type 'type'
    void registerHpxActions(); ///> Register all HPX actions

    static int main_handler(char ***argv, size_t argc);
    static int clear_handler();
    static int setNeurons_handler(const int nargs, const void *args[], const size_t []);
    static int setInputParams_handler(const Input::InputParams * data, const size_t);
    static int setMechanisms_handler (const int nargs, const void *args[], const size_t[]);
} ;
