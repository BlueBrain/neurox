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

#include "neurox/neurox.h"

static hpx_action_t main_hpx = 0;
static int main_hpx_handler( const int argc, char ** argv)
{
    //populate InputParams from command line, and broadcasts to all compute nodes
    InputParams inputParams(argc, argv);
    printf("Broadcasting InputParams...\n");
    int e = hpx_bcast_rsync(InputParams::initialize, &inputParams, sizeof (InputParams));
    assert(e == HPX_SUCCESS);

    //reads morphology data
    CoreNeuronDataLoader::loadData(argc, argv);

    nrn_finitialize( 1, inputParams.voltage );

    report_mem_usage( "After nrn_finitialize" );

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

    /// Solver execution
    BBS_netpar_solve( inputParams.tstop );

    // Report global cell statistics
    report_cell_stats();

    // prcellstate after end of solver
    //if ( globalInfo->prcellgid >= 0 ) {
    //    sprintf( prcellname, "t%g", t );
    //    prcellstate( globalInfo->prcellgid, prcellname );
    //}

    // write spike information to input_params.outpath
    output_spikes( inputParams.outputPath );

    hpx_exit(HPX_SUCCESS);
}

int main1_hpx(int argc, char** argv)
{
    //hpx initialisation
    if (hpx_init(&argc, &argv) != 0)
    {
        printf("HPX failed to initialize!\n");
        return -1;
    }
    printf("\nHPX started. Localities: %d, Threads/locality: %d\n", HPX_LOCALITIES, HPX_THREADS);

    //register HPX methods
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED, main_hpx, main_hpx_handler, HPX_INT, HPX_POINTER);
    InputParams::registerHpxActions();
    Brain::registerHpxActions();
    Neuron::registerHpxActions();
    Branch::registerHpxActions();

    //start HPX
    int e = hpx_run(&main_hpx, argc, argv);

    //clean up
    hpx_finalize();
    return e;
}
