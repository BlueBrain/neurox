#pragma once

#include "neurox/Neurox.h"

namespace Neurox
{

/**
 * @brief The Neuron class
 * Represents a neuron as its metadata and a tree-based orphology
 */
class Neuron
{
  public:
    Neuron() = delete;
    ~Neuron();

    //neuron metadata
    int id;					///> neuron global id
    hpx_t soma;		///> hpx address of the top compartment (soma)
    //TODO this should be a Branch of size 1 instead, no need for hpx_t
    
    double t; ///> current time step

    //from NrnThread
    double cj; ///<1st or 2nd order solver ... (?)
    
    //outgoing synapses
    double thresholdAP;     ///> Action Potential threshold
    int synapsesCount;      ///> number of outgoing synapses
    Synapse * synapses;       ///> outgoing Synapses (branches addr)

    static void registerHpxActions(); ///> Register all HPX actions
    static void setupTreeMatrixMinimal(Neuron * local); ///>set_tree_matrix_minimal
    static hpx_t fireActionPotential(Neuron * local); ///> fires AP, returns LCO for sent synapses
    static hpx_action_t finitialize;  ///> finitialize.c
    static hpx_action_t init;         ///> Initializes Neuron
    static hpx_action_t solve;        ///> Main loop, the solver
    static hpx_action_t addSynapses;  ///> Inserts Synapses (targets) in this Neuron

  private:
    static int finitialize_handler(); ///> initialize.c
    static int solve_handler(); ///> BBS_netpar_solve( inputParams.tstop );
    static int init_handler(const int gid, const hpx_t soma, double thresholdAP); ///> HPX constructor
    static int addSynapses_handler(Synapse * synapses, size_t synapsesCount); ///>Inserts Synapses
};

}
