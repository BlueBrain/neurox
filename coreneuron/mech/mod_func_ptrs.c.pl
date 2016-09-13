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


# Construct the modl_reg() function from a provided list
# of modules.

# Usage: mod_func.c.pl [MECH1.mod MECH2.mod ...]

@mods=@ARGV;
s/\.mod$// foreach @mods;

@mods=sort @mods;

@funcs=('init','cur','jacob','state');

print <<"__eof";
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include \"coreneuron/nrnoc/membfunc.h\"
#include \"coreneuron/nrnoc/multicore.h\"
__eof

#Get the correct SUFFIX from each mod file for each mechanism
@suffixes=();

for $m(@mods) {
  $filename = "${m}.mod";
  open $fh, '<', $filename or die "error: unable to open file '$filename' for reading : $!. ERROR1";
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
    die "error: unable to find mechanism name for ${filename}. Add the missing keywork to coreneuron/mech/mod_func.c.pl. ERROR2"
  }
  @lines[0] =~ s/^\s+|\s+$//g;  #remove trailing whitespaces from beginning and end
  @lines[0] =~ s/ +/ /; #replace multiple spaces by one
  @lines[0] =~ s/[\r\n]+$//; #remove bad endings (breakline)
  my @words = split / /, @lines[0]; #get words from first (and only) line containing 'SUFFIX'
  my $suffix = @words[1]; #get SUFFIX name as second word"
  push(@suffixes, $suffix); 
}

#Output the get of function pointers for init, jacob, current and state functions

for $f(@funcs) {
print <<"__eof";

extern void \n  @{[join ",\n  ", map {"_nrn_${f}__${_}(NrnThread*, Memb_list*, int)"} @suffixes]};
__eof
}

for $f(@funcs) {
print <<"__eof";

mod_f_t get_${f}_function(const char * sym) {

@{[join "\n",map {"  if (strcmp(sym, \"${_}\") == 0)  return _nrn_${f}__${_};"} @suffixes]}

  printf(\"Could not find ${f} function for mechanism %s\\n\", sym);
  abort();
  return NULL;
}
__eof

}

print <<"__eof";

mod_f_t get_BA_function(const char * sym, int BA_func_id)
{
  return NULL;
}
__eof

