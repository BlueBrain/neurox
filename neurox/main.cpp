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
static int main_hpx_handler( char ** argv, const int argc)
{
    //populate InputParams with info past from command line
    InputParams inputParams(argc, argv);

    //memory footprint after HPX initialisation
    report_mem_usage( "After hpx_init" );

    /////// OLD CORENEURON CODE --- WILL GO AWAY at some point ////

    char prcellname[1024], filesdat_buf[1024], datpath[1024];

    // initialise default coreneuron parameters
    //initnrn(); //part of GlobalInfo constructor

    // handles coreneuron configuration parameters
    cn_input_params input_params;

    // read command line parameters
    input_params.read_cb_opts( argc, argv );

    // set global variables for start time, timestep and temperature
    t = inputParams.tstart; // input_params.tstart;
    dt = inputParams.dt ; //input_params.dt;
    rev_dt = inputParams.rev_dt; //(int)(1./dt);
    celsius = inputParams.celsius ; //input_params.celsius;

    // full path of files.dat file
    sd_ptr filesdat=input_params.get_filesdat_path(filesdat_buf,sizeof(filesdat_buf));

    // memory footprint after mpi initialisation
    report_mem_usage( "After MPI_Init" );

    // reads mechanism information from bbcore_mech.dat
    mk_mech( datpath );

    report_mem_usage( "After mk_mech" );

    // create net_cvode instance
    mk_netcvode();

    // One part done before call to nrn_setup. Other part after.
    if ( inputParams.patternStim ) {
        nrn_set_extra_thread0_vdata();
    }

    report_mem_usage( "Before nrn_setup" );

    // reading *.dat files and setting up the data structures
    nrn_setup( input_params, filesdat, nrn_need_byteswap);

    report_mem_usage( "After nrn_setup " );

    // Invoke PatternStim
    if ( inputParams.patternStim) {
        nrn_mkPatternStim( input_params.patternstim );
    }

    /// Setting the timeout
    nrn_set_timeout(200.);

    // show all configuration parameters for current run
    input_params.show_cb_opts();

    //////////// END OF CORENEURON OLD CODE ///////////////////

    //We take all CoreNeuron data type.size() convert to hpx based data types
    //TODO: in the future we want an implementation to run from CoreNeuron
    //and one to be a stand-alone app (read directly from file)

    //Broadcast input arguments, and initalizes
    printf("Broadcasting InputParams...\n");
    int e = hpx_bcast_rsync(InputParams::initialize, &inputParams, sizeof (InputParams));
    assert(e == HPX_SUCCESS);

    CoreNeuronDataLoader::loadData();

    //Clean core neuron data, work only with HPX data
    //nrn_cleanup();

    nrn_finitialize( 1, input_params.voltage );

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
    BBS_netpar_solve( input_params.tstop );

    // Report global cell statistics
    report_cell_stats();

    // prcellstate after end of solver
    //if ( globalInfo->prcellgid >= 0 ) {
    //    sprintf( prcellname, "t%g", t );
    //    prcellstate( globalInfo->prcellgid, prcellname );
    //}

    // write spike information to input_params.outpath
    output_spikes( input_params.outpath );

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
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED, main_hpx, main_hpx_handler, HPX_POINTER, HPX_INT);
    InputParams::registerHpxActions();
    Brain::registerHpxActions();
    Neuron::registerHpxActions();
    Branch::registerHpxActions();

    //start HPX
    int e = hpx_run(&main_hpx, argv, argc);

    //clean up
    hpx_finalize();
    return e;
}
