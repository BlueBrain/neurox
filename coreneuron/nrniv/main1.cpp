/*
Copyright (c) 2014 EPFL-BBP, All rights reserved.

THIS SOFTWARE IS PROVIDED BY THE BLUE BRAIN PROJECT "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE BLUE BRAIN PROJECT
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/**
 * @file main1.cpp
 * @date 26 Oct 2014
 * @brief File containing main driver routine for CoreNeuron
 */

#include "coreneuron/nrnconf.h"
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrnoc/nrnoc_decl.h"
#include "coreneuron/nrnmpi/nrnmpi.h"
#include "coreneuron/nrniv/nrniv_decl.h"
#include "coreneuron/nrniv/output_spikes.h"
#include "coreneuron/utils/endianness.h"
#include "coreneuron/utils/memory_utils.h"
#include "coreneuron/nrniv/nrnoptarg.h"
#include "coreneuron/utils/randoms/nrnran123.h"
#include "coreneuron/utils/sdprintf.h"
#include "coreneuron/nrniv/nrn_stats.h"

#include "neurox/neurox.h"
#include "neurox/nrn_setup.h"

static hpx_action_t main_hpx = 0;
static int main_hpx_handler( cn_input_params * input_params_ptr, const size_t size )
{
    cn_input_params input_params = *input_params_ptr;

    //We take all CoreNeuron data types and convert to hpx based data types
    nrn_setup_hpx();

    //Clean core neuron data, work only with HPX data
    //nrn_cleanup();

    nrn_finitialize( 1, input_params.voltage );

    report_mem_usage( "After nrn_finitialize" );

    // call prcellstae for prcellgid
    char prcellname[1024];
    if ( input_params.prcellgid >= 0 ) {
        sprintf( prcellname, "t%g", t );
        prcellstate( input_params.prcellgid, prcellname );
    }

    // handle forwardskip
    if ( input_params.forwardskip > 0.0 ) {
        handle_forward_skip( input_params.forwardskip, input_params.prcellgid );
    }

    /// Solver execution
    BBS_netpar_solve( input_params.tstop );

    // Report global cell statistics
    report_cell_stats();

    // prcellstate after end of solver
    if ( input_params.prcellgid >= 0 ) {
        sprintf( prcellname, "t%g", t );
        prcellstate( input_params.prcellgid, prcellname );
    }

    // write spike information to input_params.outpath
    output_spikes( input_params.outpath );

    hpx_exit(HPX_SUCCESS);
}

int main1( int argc, char **argv, char **ev )
{
    //hpx initialisation
    if (hpx_init(&argc, &argv) != 0)
    {
        printf("HPX failed to initialize!\n");
        return -1;
    }
    printf("\nHPX started. Localities: %d, Threads/locality: %d\n", HPX_LOCALITIES, HPX_THREADS);

    // initialise default coreneuron parameters
    initnrn();

    // handles coreneuron configuration parameters
    cn_input_params input_params;

    // read command line parameters
    input_params.read_cb_opts( argc, argv );

    char filesdat_buf[1024];

    // set global variables for start time, timestep and temperature
    t = input_params.tstart;
    dt = input_params.dt;
    rev_dt = (int)(1./dt);
    celsius = input_params.celsius;

    // full path of files.dat file
    sd_ptr filesdat=input_params.get_filesdat_path(filesdat_buf,sizeof(filesdat_buf));

    // memory footprint after HPX initialisation
    report_mem_usage( "After hpx_init" );

    // reads mechanism information from bbcore_mech.dat
    mk_mech( input_params.datpath );

    report_mem_usage( "After mk_mech" );

    // create net_cvode instance
    mk_netcvode();

    // One part done before call to nrn_setup. Other part after.
    if ( input_params.patternstim ) {
        nrn_set_extra_thread0_vdata();
    }

    report_mem_usage( "Before nrn_setup" );

    // reading *.dat files and setting up the data structures
    //TODO: ask how to get ALL nodes to read ALL neurons (otherwise the hpx implementation doesn't work)
    nrn_setup( input_params, filesdat, nrn_need_byteswap);

    report_mem_usage( "After nrn_setup " );

    // Invoke PatternStim
    if ( input_params.patternstim ) {
        nrn_mkPatternStim( input_params.patternstim );
    }

    /// Setting the timeout
    nrn_set_timeout(200.);

    // show all configuration parameters for current run
    input_params.show_cb_opts();

    //register HPX methods
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED, main_hpx, main_hpx_handler, HPX_POINTER, HPX_SIZE_T);
    void nrn_setup_register_hpx_actions();

    //start HPX
    int e = hpx_run(&main_hpx,  &input_params, sizeof(input_params));

    //clean up
    hpx_finalize();
    return e;
}


/* perform forwardskip and call prcellstate for prcellgid */
void handle_forward_skip( double forwardskip, int prcellgid )
{
    double savedt = dt;
    double savet = t;

    dt = forwardskip * 0.1;
    t = -1e9;

    for ( int step = 0; step < 10; ++step ) {
        nrn_fixed_step_minimal();
    }

    if ( prcellgid >= 0 ) {
        prcellstate( prcellgid, "fs" );
    }

    dt = savedt;
    t = savet;
    dt2thread(-1.);
}


const char *nrn_version( int )
{
    return "version id unimplemented";
}
