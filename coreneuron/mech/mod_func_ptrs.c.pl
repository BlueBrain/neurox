#!/usr/bin/perl

# Author: Bruno Magalhaes
# Usage: mod_func_ptrs.c.pl [PATH-MECH1.mod PATH-MECH2.mod ...]

@mods=@ARGV;

s/\.mod$// foreach @mods;

@mods=sort @mods;

@funcs=('init','state');
#@funcs=('init','cur','jacob', 'state');

print <<"__eof";
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "coreneuron/coreneuron.h"
namespace coreneuron {
__eof

#Get the correct SUFFIX from each mod file for each mechanism
@suffixes_all=();
@suffixes_with_cur=(); #with cur function (BREAKPOINT block in mod)
@suffixes_with_state_vars=(); # with state vars i.e. STATE block in mod
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

  #add only those with nrn_cur function definition
  my @breakpointlines = grep /STATE/, @content;
  if (scalar @breakpointlines == 1) {
    push(@suffixes_with_state_vars, $suffix);
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
elsif ($f eq "state_vars"){
  @suffixes_with_this_func = @suffixes_with_state_vars;
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

extern void \n  @{[join ",\n  ", map {"_net_receive__${_}(NrnThread *, Memb_list*, int, int, double)"} @suffixes_with_net_receive]};

pnt_receive_t get_net_receive_function(const char * sym)
{
@{[join "\n",map {"  if (strcmp(sym, \"${_}\") == 0)  return _net_receive__${_};"} @suffixes_with_net_receive]}
  return NULL;
}

__eof

print <<"__eof";

extern void \n  @{[join ",\n  ", map {"_nrn_cur__${_}(NrnThread*, Memb_list*, int, mod_acc_f_t = NULL, mod_acc_f_t = NULL, void * args = NULL)"} @suffixes_with_cur]};

mod_cur_f_t get_cur_function(const char * sym)
{
@{[join "\n",map {"  if (strcmp(sym, \"${_}\") == 0)  return _nrn_cur__${_};"} @suffixes_with_cur]}
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

#output jacob functions (TODO not available yet)
print <<"__eof";

mod_f_t get_jacob_function(const char * sym)
{
  return NULL;
}

__eof


# CVODES-specifc methods
print <<"__eof";

extern void \n  @{[join ",\n  ", map {"_nrn_ode_state_vars__${_}(int*, int**, int**)"} @suffixes_with_state_vars]};

state_vars_f_t get_ode_state_vars_function(const char * sym)
{
@{[join "\n",map {"  if (strcmp(sym, \"${_}\") == 0)  return _nrn_ode_state_vars__${_};"} @suffixes_with_state_vars]}
  return NULL;
}


extern int \n  @{[join ",\n  ", map {"_nrn_ode_matsol1__${_}(threadargsproto)"} @suffixes_with_state_vars]};

cvode_f_t get_ode_matsol_function(const char * sym)
{
@{[join "\n",map {"  if (strcmp(sym, \"${_}\") == 0)  return _nrn_ode_matsol1__${_};"} @suffixes_with_state_vars]}
  return NULL;
}


extern int \n  @{[join ",\n  ", map {"_nrn_ode_spec1__${_}(threadargsproto)"} @suffixes_with_state_vars]};

cvode_f_t get_ode_spec_function(const char * sym)
{
@{[join "\n",map {"  if (strcmp(sym, \"${_}\") == 0)  return _nrn_ode_spec1__${_};"} @suffixes_with_state_vars]}
  return NULL;
}

} //namespace coreneuron
__eof

