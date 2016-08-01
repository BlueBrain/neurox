#include "neurox/Neurox.h"
#include <cstring>

#include "nrniv/nrniv_decl.h"
#include "nrniv/nrn_stats.h"

namespace Neurox
{

int neuronsCount=-1;
hpx_t neuronsAddr = HPX_NULL;
int mechanismsCount=-1;
extern int * mechanismsMap = nullptr;
Neurox::Mechanism * mechanisms = nullptr;
Input::InputParams * inputParams = nullptr;

hpx_action_t setNeurons = 0;
int setNeurons_handler(const int nargs, const void *args[], const size_t[])
{
    /** nargs=2 where
     * args[0] = neuronsCount
     * args[1] = neuronsAddr
     */

    neurox_hpx_pin(uint64_t);
    assert(nargs==2);
    neuronsCount = *(int*)args[0];
    neuronsAddr = *(hpx_t*)args[1];
    neurox_hpx_unpin;
}

hpx_action_t setInputParams = 0;
int setInputParams_handler(const Input::InputParams * data, const size_t)
{
    neurox_hpx_pin(uint64_t);
    if (inputParams!=nullptr)
        delete [] Neurox::inputParams;

    inputParams = new Neurox::Input::InputParams();
    memcpy(Neurox::inputParams, data, sizeof(Neurox::Input::InputParams));
    neurox_hpx_unpin;
}

hpx_action_t setMechanisms = 0;
int setMechanisms_handler(const int nargs, const void *args[], const size_t sizes[])
{
    /**
     * nargs=3 where:
     * args[0] = array of all mechanisms info
     * args[1] = array of all mechanisms dependencies (children in mechanisms tree)
     * args[2] = array of all mechanisms names (sym)
     */

    neurox_hpx_pin(uint64_t);
    assert(nargs==3);
    mechanismsCount = sizes[0]/sizeof(Mechanism);
    mechanisms = new Mechanism[mechanismsCount];

    int offsetChildren=0;
    int offsetSym=0;
    int maxMechType=-1;
    for (int m=0; m<mechanismsCount; m++)
    {
        Mechanism & mech = ((Mechanism*) args[0])[m];
        int * children = mech.childrenCount == 0 ? nullptr : &((int*) args[1])[offsetChildren];
        char * sym = mech.symLength == 0 ? nullptr : &((char*) args[2])[offsetSym];
        Mechanism * finalMech = new Mechanism(
                    mech.type, mech.dataSize, mech.pdataSize,
                    mech.isArtificial, mech.pntMap, mech.isIon,
                    mech.symLength, sym,
                    mech.isTopMechanism, mech.childrenCount, children);
        mechanisms[m] = *finalMech;
        offsetChildren +=  mech.childrenCount;
        offsetSym += mech.symLength;
        if (mech.type > maxMechType)
            maxMechType = mech.type;
    }

    //initializes map of mechanisms ids to offset
    mechanismsMap = new int[maxMechType];
    for (int i=0; i<maxMechType; i++)
        mechanismsMap[i]=-1;
    for (int m=0; m<mechanismsCount; m++)
        mechanismsMap[mechanisms[m].type]=m;

    neurox_hpx_unpin;
}

hpx_action_t main = 0;
static int main_handler( char **argv, size_t argc)
{
    //reads morphology data
    Neurox::Input::Coreneuron::DataLoader::loadData(argc, argv);

    //populate InputParams from command line, and broadcasts to all compute nodes
    Input::InputParams inputParams(argc, argv);
    printf("Broadcasting InputParams...\n");
    int e = hpx_bcast_rsync(Neurox::setInputParams, &inputParams, sizeof (Input::InputParams));
    assert(e == HPX_SUCCESS);

    //call finitialize.c (nrn_finitialize( 1, inputParams.voltage )
    hpx_par_for_sync( [&] (int i, void*) -> int
        { return hpx_call_sync(getNeuronAddr(i), Neuron::finitialize, NULL, 0);},
        0, neuronsCount, NULL);

    // call prcellstae for prcellgid
    //opens the file that will store this cell's info
    //if ( globalInfo->prcellgid >= 0 ) {
    //    sprintf( prcellname, "t%g", t );
    //    prcellstate( globalInfo->prcellgid, prcellname );
    //}

    // handle forwardskip
    // many steps with large dt so that cells start at their resting potential
    //if ( input_params.forwardskip > 0.0 ) {
    //    handle_forward_skip( input_params.forwardskip, input_params.prcellgid );
    //}

    Solver::BackwardEuler::solve(inputParams.dt, inputParams.tstop); //BBS_netpar_solve( inputParams.tstop );
    assert(e == HPX_SUCCESS);

    BBS_netpar_solve( inputParams.tstop );

    // Report global cell statistics
    report_cell_stats();

    // prcellstate after end of solver
    //if ( globalInfo->prcellgid >= 0 ) {
    //    sprintf( prcellname, "t%g", t );
    //    prcellstate( globalInfo->prcellgid, prcellname );
    //}

    // write spike information to input_params.outpath
    //output_spikes( inputParams.outputPath );

    hpx_exit(HPX_SUCCESS);
}

void registerHpxActions()
{
    neurox_hpx_register_action(1,Neurox::main);
    neurox_hpx_register_action(2,Neurox::setNeurons);
    neurox_hpx_register_action(1,Neurox::setInputParams);
    neurox_hpx_register_action(2,Neurox::setMechanisms);

}

};
