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

#include <stdlib.h>
#include <math.h>
#include "coreneuron/nrnconf.h"
#include "coreneuron/utils/randoms/nrnran123.h"
#include "coreneuron/utils/randoms/Random123/philox.h"
#include "coreneuron/nrniv/nrnmutdec.h"

static const double SHIFT32 = 1.0 / 4294967297.0; /* 1/(2^32 + 1) */

static philox4x32_key_t k = {{0}};

static size_t instance_count_ = 0;
size_t nrnran123_instance_count() {
    return instance_count_;
}

/* now this is declated in nrnran123.h so that its available for prototype declaration

struct nrnran123_State {
        philox4x32_ctr_t c;
        philox4x32_ctr_t r;
        unsigned char which_;
};
*/

size_t nrnran123_state_size() {
    return sizeof(nrnran123_State);
}

void nrnran123_set_globalindex(uint32_t gix) {
    k.v[0] = gix;
}

/* if one sets the global, one should reset all the stream sequences. */
uint32_t nrnran123_get_globalindex() {
    return k.v[0];
}

static MUTDEC void nrnran123_mutconstruct() {
    if (!mut_) {
        MUTCONSTRUCT(1);
    }
}

nrnran123_State* nrnran123_newstream(uint32_t id1, uint32_t id2) {
    return nrnran123_newstream3(id1, id2, 0);
}

nrnran123_State* nrnran123_newstream3(uint32_t id1, uint32_t id2, uint32_t id3) {
    nrnran123_State* s = (nrnran123_State*)ecalloc(sizeof(nrnran123_State), 1);
    s->c.v[1] = id3;
    s->c.v[2] = id1;
    s->c.v[3] = id2;
    nrnran123_setseq(s, 0, 0);
    MUTLOCK
    ++instance_count_;
    MUTUNLOCK
    return s;
}

void nrnran123_deletestream(nrnran123_State* s) {
    MUTLOCK
    --instance_count_;
    MUTUNLOCK
    free(s);
}

void nrnran123_getseq(nrnran123_State* s, uint32_t* seq, unsigned char* which) {
    *seq = s->c.v[0];
    *which = s->which_;
}

void nrnran123_setseq(nrnran123_State* s, uint32_t seq, unsigned char which) {
    if (which > 3) {
        s->which_ = 0;
    } else {
        s->which_ = which;
    }
    s->c.v[0] = seq;
    s->r = philox4x32(s->c, k);
}

void nrnran123_getids(nrnran123_State* s, uint32_t* id1, uint32_t* id2) {
    *id1 = s->c.v[2];
    *id2 = s->c.v[3];
}

uint32_t nrnran123_ipick(nrnran123_State* s) {
    uint32_t rval;
    unsigned char which = s->which_;
    assert(which < 4);
    rval = s->r.v[which++];
    if (which > 3) {
        which = 0;
        s->c.v[0]++;
        s->r = philox4x32(s->c, k);
    }
    s->which_ = which;
    return rval;
}

double nrnran123_dblpick(nrnran123_State* s) {
    return nrnran123_uint2dbl(nrnran123_ipick(s));
}

double nrnran123_negexp(nrnran123_State* s) {
    /* min 2.3283064e-10 to max 22.18071 */
    return -log(nrnran123_dblpick(s));
}

/* at cost of a cached  value we could compute two at a time. */
double nrnran123_normal(nrnran123_State* s) {
    double w, x, y;
    double u1, u2;
    do {
        u1 = nrnran123_dblpick(s);
        u2 = nrnran123_dblpick(s);
        u1 = 2. * u1 - 1.;
        u2 = 2. * u2 - 1.;
        w = (u1 * u1) + (u2 * u2);
    } while (w > 1);

    y = sqrt((-2. * log(w)) / w);
    x = u1 * y;
    return x;
}

double nrnran123_uint2dbl(uint32_t u) {
    /* 0 to 2^32-1 transforms to double value in open (0,1) interval */
    /* min 2.3283064e-10 to max (1 - 2.3283064e-10) */
    return ((double)u + 1.0) * SHIFT32;
}
