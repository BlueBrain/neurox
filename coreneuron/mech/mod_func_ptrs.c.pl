#!/usr/bin/perl

# Copyright (c) 2014 EPFL-BBP, All rights reserved.
# 
# THIS SOFTWARE IS PROVIDED BY THE BLUE BRAIN PROJECT "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE BLUE BRAIN PROJECT
# BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
# IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


# Constructs the getters for mechanisms function pointers;

# Usage: mod_func_ptrs.c.pl [PATH-MECH1.mod PATH-MECH2.mod ...]

@mods=@ARGV;

s/\.mod$// foreach @mods;

@mods=sort @mods;

@funcs=('init','cur','state');
#@funcs=('init','cur','jacob', 'state');

print <<"__eof";
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "coreneuron/coreneuron.h"
__eof

#Get the correct SUFFIX from each mod file for each mechanism
@suffixes_all=();
@suffixes_with_cur=(); #with cur function (BREAKPOINT block in mod)
@suffixes_with_net_receive=(); #with NET_RECEIVE function in mod

for $m(@mods) {
  $filename = "${m}.mod";
  open $fh, '<', $filename or die "error: unable to open file '$filename' for reading : $!.";
  my @content = <$fh>;
  close $fh;
  my @lines = grep /SUFFIX/, @content;
  if (!@lines) {
    @lines = grep /POINT_PROCESS/, @content;
  }
  if (!@lines) {
    @lines = grep /ARTIFICIAL_CELL/, @content;
  }
  if (!@lines) {
    die "error: unable to find mechanism name for ${filename}. Add the missing keywork to coreneuron/mech/mod_func.c.pl."
  }
  @lines[0] =~ s/^\s+|\s+$//g;  #remove trailing whitespaces from beginning and end
  @lines[0] =~ s/ +/ /; #replace multiple spaces by one
  @lines[0] =~ s/[\r\n]+$//; #remove bad endings (breakline)
  my @words = split / /, @lines[0]; #get words from first (and only) line containing 'SUFFIX'
  my $suffix = @words[1]; #get SUFFIX name as second word"
  push(@suffixes_all, $suffix);

  #add only those with nrn_cur function definition
  my @breakpointlines = grep /BREAKPOINT/, @content;
  if (scalar @breakpointlines == 1) {
    push(@suffixes_with_cur, $suffix);
  }

  #add only those with net_receive function definition
  my @breakpointlines = grep /NET_RECEIVE/, @content;
  if (scalar @breakpointlines == 1) {
    push(@suffixes_with_net_receive, $suffix);
  }
}

#Output the get of function pointers for init, jacob, current and  state functions

for $f(@funcs) {

@suffixes_with_this_func=();
if ($f eq "cur"){
  @suffixes_with_this_func = @suffixes_with_cur;
}
else {
  @suffixes_with_this_func = @suffixes_all;
}

print <<"__eof";

extern void \n  @{[join ",\n  ", map {"_nrn_${f}__${_}(NrnThread*, Memb_list*, int)"} @suffixes_with_this_func]};

mod_f_t get_${f}_function(const char * sym)
{
@{[join "\n",map {"  if (strcmp(sym, \"${_}\") == 0)  return _nrn_${f}__${_};"} @suffixes_with_this_func]}
  return NULL;
}

__eof
}

#Output the get of function pointers for init, jacob, current, state and net_receive functions

print <<"__eof";

extern void \n  @{[join ",\n  ", map {"_net_receive2__${_}(Point_process*, int, double)"} @suffixes_with_net_receive]};

pnt_receive2_t get_net_receive_function(const char * sym)
{
@{[join "\n",map {"  if (strcmp(sym, \"${_}\") == 0)  return _net_receive2__${_};"} @suffixes_with_net_receive]}
  return NULL;
}

__eof

# output of locked current function
#ions only: 'Nap_Et2','NaTa_t','NaTs2_t','SKv3_1','Im','K_Pst', 'K_Tst', 'SK_E2', 'Ca', 'CaDynamics_E2', 'Ca_HVA', 'Ca_LVAst';

print <<"__eof";

extern void \n  @{[join ",\n  ", map {"_nrn_cur_parallel__${_}(NrnThread*, Memb_list*, int)"} @suffixes_with_cur]};

mod_parallel_f_t get_cur_parallel_function(const char * sym)
{
@{[join "\n",map {"  if (strcmp(sym, \"${_}\") == 0)  return _nrn_cur_parallel__${_};"} @suffixes_with_cur]}
  return NULL;
}
__eof


#output get_state_vars functions (ions only)

print <<"__eof";

extern void \n  @{[join ",\n  ", map {"_nrn_stave_vars__${_}(short*, short**, short**)"} @suffixes_with_cur]};

state_vars_f_t get_state_vars_function(const char * sym)
{
@{[join "\n",map {"  if (strcmp(sym, \"${_}\") == 0)  return _nrn_stave_vars__${_};"} @suffixes_with_cur]}
  return NULL;
}
__eof

#output BA functions (not available yet)
print <<"__eof";

mod_f_t get_BA_function(const char * sym, int BA_func_id)
{
  return NULL;
}
__eof

#output jacob functions (not available yet)
# TODO get rid of this
print <<"__eof";
mod_f_t get_jacob_function(const char * sym)
{
  return NULL;
}
__eof


