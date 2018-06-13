#!/usr/bin/python

# read the translated mod files and construct _kinderiv.h
# this is implemented to avoid function callbacks using
# function pointers while using methods like euler, newton
# or kinetic schemes. This is required for GPU implementation.

# if _kinderiv.h already exists and is the same as what this script generates
# then do not update it (to avoid re-compiling all the mech c files if
# only one of them changes)

kftmp = '_tmp_kinderiv.h'
kf = '_kinderiv.h'

import os

fnames = [f.replace('.mod', '.cpp') for f in os.listdir('.') if f.endswith('.mod')]
euler = []
deriv = []
kin = []

for fname in fnames:
  f = open(fname, "r")
  for line in f:
    word = line.split()
    if len(word) > 3:
      if word[0] == '/*' and word[1] == '_derivimplicit_':
        deriv.append([word[2], word[3], fname, word[1]])
      if word[0] == '/*' and word[1] == '_kinetic_':
        kin.append([word[2], word[3], fname, word[1]])
      if word[0] == '/*' and word[1] == '_euler_':
        euler.append([word[2], word[3], fname, word[1]])
  f.close()

fout = open(kftmp, "w")
fout.write('''
#ifndef _kinderiv_h
#define _kinderiv_h
''')

fout.write("\n/* data used to construct this file */\n")
for l in [deriv, kin,  euler]:
  for item in l:
    fout.write('/*')
    for word in item:
      fout.write(' %s' % word)
    fout.write(' */\n')

fout.write("\n/* declarations */\n")
fout.write("\nnamespace coreneuron {\n")
for item in deriv:
  fout.write('#pragma acc routine seq\n')
  fout.write('extern int %s%s(_threadargsproto_);\n' % (item[0], item[1]))
  fout.write('#pragma acc routine seq\n')
  fout.write('extern int _newton_%s%s(_threadargsproto_);\n' % (item[0], item[1]))

for item in kin:
  fout.write('#pragma acc routine seq\n')
  fout.write('extern int %s%s(void*, double*, _threadargsproto_);\n' % (item[0], item[1]))

for item in euler:
  fout.write('#pragma acc routine seq\n')
  fout.write('extern int %s%s(_threadargsproto_);\n' % (item[0], item[1]))

fout.write("\n/* callback indices */\n")
derivoffset = 1
kinoffset = 1

for i, item in enumerate(deriv):
  fout.write('#define _derivimplicit_%s%s %d\n' % (item[0], item[1], i + derivoffset))
for i, item in enumerate(kin):
  fout.write('#define _kinetic_%s%s %d\n' % (item[0], item[1], i + kinoffset))
for i, item in enumerate(euler):
  fout.write('#define _euler_%s%s %d\n' % (item[0], item[1], i + derivoffset))

fout.write("\n/* switch cases */\n")
fout.write("\n#define _NRN_DERIVIMPLICIT_CASES \\\n")
for item in (deriv):
  fout.write("  case _derivimplicit_%s%s: %s%s(_threadargs_); break; \\\n" % (item[0], item[1], item[0], item[1]))
fout.write("\n")

fout.write("\n#define _NRN_DERIVIMPLICIT_NEWTON_CASES \\\n")
for item in deriv:
  fout.write("  case _derivimplicit_%s%s: _newton_%s%s(_threadargs_); break; \\\n" % (item[0], item[1], item[0], item[1]))
fout.write("\n")

fout.write("\n#define _NRN_KINETIC_CASES \\\n")
for item in kin:
  fout.write("  case _kinetic_%s%s: %s%s(so, rhs, _threadargs_); break; \\\n" % (item[0], item[1], item[0], item[1]))
fout.write("\n")

fout.write("\n#define _NRN_EULER_CASES \\\n")
for item in euler:
  fout.write("  case _euler_%s%s: %s%s(_threadargs_); break; \\\n" % (item[0], item[1], item[0], item[1]))
fout.write("\n")

fout.write('\n#endif\n')

fout.write("\n} //namespace coreneuron\n")
fout.close()

# if kf exists and is same as kftmp, just remove kftmp. Otherwise
# rename kftmp to kf
import filecmp
b = False
try:
  b = filecmp.cmp(kftmp, kf)
except:
  pass
if b:
  os.remove(kftmp)
else:
  os.rename(kftmp, kf)
