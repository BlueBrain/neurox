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

#include <iostream>
#include "coreneuron/nrnconf.h"
#include "coreneuron/nrniv/nrniv_decl.h"
#include "coreneuron/nrniv/output_spikes.h"
#include "coreneuron/nrnmpi/nrnmpi.h"
#include "coreneuron/nrniv/nrnmutdec.h"
#include "coreneuron/utils/sdprintf.h"

int spikevec_buffer_size;
int spikevec_size;
double* spikevec_time;
int* spikevec_gid;

static MUTDEC

    void
    mk_spikevec_buffer(int sz) {
    spikevec_buffer_size = sz;
    spikevec_size = 0;
    spikevec_time = new double[sz];
    spikevec_gid = new int[sz];
    MUTCONSTRUCT(1);
}

void spikevec_lock() {
    MUTLOCK
}
void spikevec_unlock() {
    MUTUNLOCK
}

void output_spikes(const char* outpath) {
    char fnamebuf[100];
    sd_ptr fname = sdprintf(fnamebuf, sizeof(fnamebuf), "%s/out%d.dat", outpath, nrnmpi_myid);
    FILE* f = fopen(fname, "w");
    if (!f && nrnmpi_myid == 0) {
        std::cout << "WARNING: Could not open file for writing spikes." << std::endl;
        return;
    }

    for (int i = 0; i < spikevec_size; ++i)
        if (spikevec_gid[i] > -1)
            fprintf(f, "%.8g\t%d\n", spikevec_time[i], spikevec_gid[i]);

    fclose(f);
}

void validation(std::vector<std::pair<double, int> >& res) {
    for (int i = 0; i < spikevec_size; ++i)
        if (spikevec_gid[i] > -1)
            res.push_back(std::make_pair(spikevec_time[i], spikevec_gid[i]));
}
