#include "coreneuron/nrnconf.h"
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrnoc/nrnoc_decl.h"
#include "coreneuron/nrnmpi/nrnmpi.h"
#include "coreneuron/nrniv/nrniv_decl.h"
#include "coreneuron/nrniv/nrn_setup.h"
#include "coreneuron/nrniv/nrnoptarg.h"
#include "coreneuron/nrniv/output_spikes.h"
#include "coreneuron/utils/endianness.h"
#include "coreneuron/utils/memory_utils.h"

#include "coreneuron/utils/randoms/nrnran123.h"
#include "coreneuron/utils/sdprintf.h"
#include "coreneuron/nrniv/nrn_stats.h"

#include "neurox/Neurox.h"

using namespace Neurox;

static hpx_action_t main_hpx = 0;
static int main_hpx_handler( char **argv, size_t argc)
{
    //populate InputParams from command line, and broadcasts to all compute nodes
    Input::InputParams inputParams(argc, argv);
    printf("Broadcasting InputParams...\n");
    int e = hpx_bcast_rsync(Neurox::setInputParams, &inputParams, sizeof (Input::InputParams));
    assert(e == HPX_SUCCESS);

    //reads morphology data
    Neurox::Input::Coreneuron::DataLoader::loadData(argc, argv);

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

int main1_hpx(int argc, char** argv)
{ 
    //register HPX methods
    neurox_hpx_register_action(1, main_hpx);
    Neurox::registerHpxActions();
    Neuron::registerHpxActions();
    Branch::registerHpxActions();
    Solver::BackwardEuler::registerHpxActions();

    //Init HPX
    if (hpx_init(&argc, &argv) != 0)
    {
        printf("HPX failed to initialize!\n");
        return 1;
    }
    printf("\nHPX started. Localities: %d, Threads/locality: %d\n", HPX_LOCALITIES, HPX_THREADS);

    //Run main
    int e = hpx_run(&main_hpx, argv, argc);

    //clean up
    hpx_finalize();
    return e;
}
