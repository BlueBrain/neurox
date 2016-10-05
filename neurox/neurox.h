#pragma once

//Coreneuron core datatypes
#include "coreneuron/coreneuron.h"  //definitions

//Core datatypes
#include "neurox/neurox_hpx.h"
#include "neurox/input/InputParams.h"
#include "neurox/datatypes/Event.h"
#include "neurox/datatypes/VecPlayContinuous.h"
#include "neurox/datatypes/NetCon.h"
#include "neurox/datatypes/Mechanism.h"
#include "neurox/datatypes/Branch.h"
#include "neurox/datatypes/Neuron.h"

//Fixed-step Backward-Euler solver
#include "neurox/solver/HinesSolver.h"

//CoreNeuron data loader
#include "neurox/input/coreneuron/Compartment.h"
#include "neurox/input/coreneuron/DataLoader.h"

//Miscellaneous
#include "neurox/misc/Statistics.h"

//#define NDEBUG //(if active, disables assertions and debug output)

#define IClamp 7
#define ProbAMPANMDA_EMS 137
#define ProbGABAAB_EMS 139

///neurox namespace: contains global information that is copied to all localities
namespace neurox
{
    extern int neuronsCount; 	///> total neurons count in the system
    extern hpx_t neuronsAddr; 	///> hpx address of the first position of the neurons array

    extern int mechanismsCount; ///> number of mechanisms
    extern neurox::Mechanism ** mechanisms; ///> array to all existing mechanisms
    extern int * mechanismsMap; ///>map of mechanisms offset in 'mechanisms' by 'mechanism type'

    extern Input::InputParams *  inputParams; ///> Parameters parsed from command line

    extern hpx_t timeMachine;      ///> array of AND gates to control progress of individual neurons

    extern hpx_action_t main;            ///> execution starting point (called via hpx_run)
    extern hpx_action_t clear;           ///> clears all memory utilised including neurons, branches and mechanisms information
    extern hpx_action_t setNeurons;      ///> Initialized neurons and neuronsAddr global vars
    extern hpx_action_t setInputParams;  ///> Initializes InputParams
    extern hpx_action_t setMechanisms;	 ///> Initializes Mechanisms
    extern hpx_action_t setTimeMachine;     ///> Initializes timeLCOs

    hpx_t getNeuronAddr(int i);      ///> get HPX address of i-th neuron
    Mechanism * getMechanismFromType(int type); ///> returns mechanisms of type 'type'
    void registerHpxActions(); ///> Register all HPX actions

    static int main_handler(char **argv, size_t argc);
    static int clear_handler();
    static int setNeurons_handler(const int nargs, const void *args[], const size_t []);
    static int setInputParams_handler(const Input::InputParams * data, const size_t);
    static int setMechanisms_handler (const int nargs, const void *args[], const size_t[]);
    static int setTimeMachine_handler(hpx_t * timeLCOs_ptr, size_t);
} ;
