/*
Copyright (c) 2016, Blue Brain Project
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.
*/

/**
 * @file main1.cpp
 * @date 26 Oct 2014
 * @brief File containing main driver routine for CoreNeuron
 */

#include "coreneuron/utils/randoms/nrnran123.h"
#include "coreneuron/nrnconf.h"
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrnoc/nrnoc_decl.h"
#include "coreneuron/nrnmpi/nrnmpi.h"
#include "coreneuron/nrniv/nrniv_decl.h"
#include "coreneuron/nrniv/output_spikes.h"
#include "coreneuron/utils/endianness.h"
#include "coreneuron/utils/memory_utils.h"
#include "coreneuron/nrniv/nrnoptarg.h"
#include "coreneuron/utils/sdprintf.h"
#include "coreneuron/nrniv/nrn_stats.h"
#include "coreneuron/utils/reports/nrnreport.h"
#include "coreneuron/nrniv/nrn_acc_manager.h"
#include "coreneuron/nrniv/profiler_interface.h"
#include "coreneuron/nrniv/partrans.h"
#include "coreneuron/nrniv/multisend.h"
#include <string.h>

#if 0
#include <fenv.h>
#define NRN_FEEXCEPT (FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW)
int nrn_feenableexcept() {
  int result = -1;
  result = feenableexcept(NRN_FEEXCEPT);
  return result;
}
#endif

int main1(int argc, char* argv[], char** env);
void call_prcellstate_for_prcellgid(int prcellgid, int compute_gpu, int is_init);
void nrn_init_and_load_data(int argc, char* argv[], bool run_setup_cleanup=true) {
#if defined(NRN_FEEXCEPT)
    nrn_feenableexcept();
#endif

#ifdef ENABLE_SELECTIVE_PROFILING
    stop_profile();
#endif

    // mpi initialisation
#if NRNMPI
    nrnmpi_init(1, &argc, &argv);
#endif

    // memory footprint after mpi initialisation
    report_mem_usage("After MPI_Init");

    // initialise default coreneuron parameters
    initnrn();

    // create mutex for nrn123, protect instance_count_
    nrnran123_mutconstruct();

    // read command line parameters and parameter config files
    nrnopt_parse(argc, (const char**)argv);

    // set global variables
    // precedence is: set by user, globals.dat, 34.0
    celsius = nrnopt_get_dbl("--celsius");
    t = celsius; // later will read globals.dat and compare with this.

#if _OPENACC
    if (!nrnopt_get_flag("--gpu") && nrnopt_get_int("--cell-permute") == 2) {
        fprintf(
            stderr,
            "compiled with _OPENACC does not allow the combination of --cell-permute=2 and missing --gpu\n");
        exit(1);
    }
#endif

    // if multi-threading enabled, make sure mpi library supports it
#if NRNMPI
    if (nrnopt_get_flag("--threading")) {
        nrnmpi_check_threading_support();
    }
#endif

    // full path of files.dat file
    std::string filesdat(nrnopt_get_str("--datpath") + "/" + nrnopt_get_str("--filesdat"));

    // reads mechanism information from bbcore_mech.dat
    mk_mech(nrnopt_get_str("--datpath").c_str());

    // read the global variable names and set their values from globals.dat
    set_globals(nrnopt_get_str("--datpath").c_str());

    report_mem_usage("After mk_mech");

    // set global variables for start time, timestep and temperature
    t = nrnopt_get_dbl("--tstart");

    if (nrnopt_get_dbl("--dt") != -1000.) {  // command line arg highest precedence
        dt = nrnopt_get_dbl("--dt");
    } else if (dt == -1000.) {  // not on command line and no dt in globals.dat
        dt = 0.025;             // lowest precedence
    }
    nrnopt_modify_dbl("--dt", dt);

    rev_dt = (int)(1. / dt);

    if (nrnopt_get_dbl("--celsius") != -1000.) {  // command line arg highest precedence
        celsius = nrnopt_get_dbl("--celsius");
    } else if (celsius == -1000.) {  // not on command line and no celsius in globals.dat
        celsius = 34.0;              // lowest precedence
    }
    nrnopt_modify_dbl("--celsius", celsius);

    // create net_cvode instance
    mk_netcvode();

    // One part done before call to nrn_setup. Other part after.
    if (nrnopt_get_flag("--pattern")) {
        nrn_set_extra_thread0_vdata();
    }

    report_mem_usage("Before nrn_setup");

    // set if need to interleave cells
    use_interleave_permute = nrnopt_get_int("--cell-permute");
    cellorder_nwarp = nrnopt_get_int("--nwarp");
    use_solve_interleave = nrnopt_get_int("--cell-permute");

    // pass by flag so existing tests do not need a changed nrn_setup prototype.
    nrn_setup_multiple = nrnopt_get_int("--multiple");
    nrn_setup_extracon = nrnopt_get_int("--extracon");

    // multisend options
    use_multisend_ = nrnopt_get_flag("--multisend") ? 1 : 0;
    n_multisend_interval = nrnopt_get_int("--ms-subintervals");
    use_phase2_ = (nrnopt_get_int("--ms-phases") == 2) ? 1 : 0;

    // reading *.dat files and setting up the data structures, setting mindelay
    nrn_setup(filesdat.c_str(), nrn_need_byteswap, run_setup_cleanup);

    // Allgather spike compression and  bin queuing.
    nrn_use_bin_queue_ = nrnopt_get_flag("--binqueue");
    int spkcompress = nrnopt_get_int("--spkcompress");
    nrnmpi_spike_compress(spkcompress, (spkcompress ? true : false), use_multisend_);

    report_mem_usage("After nrn_setup ");

    // Invoke PatternStim
    if (nrnopt_get_flag("--pattern")) {
        nrn_mkPatternStim(nrnopt_get_str("--pattern").c_str());
    }

    /// Setting the timeout
    nrn_set_timeout(200.);

    // show all configuration parameters for current run
    nrnopt_show();

    // allocate buffer for mpi communication
    mk_spikevec_buffer(nrnopt_get_int("--spikebuf"));

    report_mem_usage("After mk_spikevec_buffer");

    if (nrnopt_get_flag("-gpu")) {
        setup_nrnthreads_on_device(nrn_threads, nrn_nthread);
    }

    if (nrn_have_gaps) {
        nrn_partrans::gap_update_indices();
    }

    // call prcellstate for prcellgid
    call_prcellstate_for_prcellgid(nrnopt_get_int("--prcellgid"), nrnopt_get_flag("-gpu"), 1);
}

void call_prcellstate_for_prcellgid(int prcellgid, int compute_gpu, int is_init) {
    char prcellname[1024];
#ifdef ENABLE_CUDA
    const char* prprefix = "cu";
#else
    const char* prprefix = "acc";
#endif

    if (prcellgid >= 0) {
        if (compute_gpu) {
            if (is_init)
                sprintf(prcellname, "%s_gpu_init", prprefix);
            else
                sprintf(prcellname, "%s_gpu_t%g", prprefix, t);
        } else {
            if (is_init)
                strcpy(prcellname, "cpu_init");
            else
                sprintf(prcellname, "cpu_t%g", t);
        }
        update_nrnthreads_on_host(nrn_threads, nrn_nthread);
        prcellstate(prcellgid, prcellname);
    }
}

int main1(int argc, char** argv, char** env) {
    (void)env; /* unused */

    // initializationa and loading functions moved to separate
    nrn_init_and_load_data(argc, argv);
    // nrnopt_get... still available until call nrnopt_delete()

    bool compute_gpu = nrnopt_get_flag("-gpu");
    #pragma acc data copyin(celsius, secondorder) if (compute_gpu)
    {
        double v = nrnopt_get_dbl("--voltage");
        nrn_finitialize(v != 1000., v);

        report_mem_usage("After nrn_finitialize");

#ifdef ENABLE_REPORTING
        ReportGenerator* r = NULL;
#endif

        // if reports are enabled using ReportingLib
        if (nrnopt_get_flag("--report")) {
#ifdef ENABLE_REPORTING
            if (nrnopt_get_int("--multiple") > 1) {
                if (nrnmpi_myid == 0)
                    printf(
                        "\n WARNING! : Can't enable reports with model duplications feature! \n");
            } else {
                r = new ReportGenerator(nrnopt_get_int("--report"), nrnopt_get_int("--tstart"),
                                        nrnopt_get_dbl("--tstop"), nrnopt_get_int("--dt"), nrnopt_get_dbl("--mindelay"),
                                        nrnopt_get_dbl("--dt_report"), nrnopt_get_str("--outpath"));
                r->register_report();
            }
#else
            if (nrnmpi_myid == 0)
                printf("\n WARNING! : Can't enable reports, recompile with ReportingLib! \n");
#endif
        }

        // call prcellstate for prcellgid
        call_prcellstate_for_prcellgid(nrnopt_get_int("--prcellgid"), compute_gpu, 0);

        // handle forwardskip
        if (nrnopt_get_dbl("--forwardskip") > 0.0) {
            handle_forward_skip(nrnopt_get_dbl("--forwardskip"), nrnopt_get_int("--prcellgid"));
        }

#ifdef ENABLE_SELECTIVE_PROFILING
        start_profile();
#endif

        /// Solver execution
        BBS_netpar_solve(nrnopt_get_dbl("--tstop"));

        // Report global cell statistics
        report_cell_stats();

#ifdef ENABLE_SELECTIVE_PROFILING
        stop_profile();
#endif

        // prcellstate after end of solver
        call_prcellstate_for_prcellgid(nrnopt_get_int("--prcellgid"), compute_gpu, 0);

#ifdef ENABLE_REPORTING
        if (nrnopt_get_int("--report") && r)
            delete r;
#endif
    }

    // write spike information to outpath
    output_spikes(nrnopt_get_str("--outpath").c_str());

    // Cleaning the memory
    nrn_cleanup();

    // mpi finalize
#if NRNMPI
    nrnmpi_finalize();
#endif

    finalize_data_on_device();

    return 0;
}

/* perform forwardskip and call prcellstate for prcellgid */
void handle_forward_skip(double forwardskip, int prcellgid) {
    double savedt = dt;
    double savet = t;

    dt = forwardskip * 0.1;
    t = -1e9;

    for (int step = 0; step < 10; ++step) {
        nrn_fixed_step_minimal();
    }

    if (prcellgid >= 0) {
        prcellstate(prcellgid, "fs");
    }

    dt = savedt;
    t = savet;
    dt2thread(-1.);
}

const char* nrn_version(int) {
    return "version id unimplemented";
}
