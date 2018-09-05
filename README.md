# neurox

A scalable parallel & distributed asynchronous simulator of extended Hodgkin-Huxley neuron models

## Installation

Dependencies
- cmake 2.8.12+: www.cmake.org
- hpx 4+: https://hpx.crest.iu.edu/
- tclap 1.2.1+: http://tclap.sourceforge.net/
- libCoreNeuron: https://github.com/BlueBrain/coreneuron
- sundials CVODES 3.1.0+: https://computation.llnl.gov/projects/sundials/cvodes
  - requires SUPERLUMT solver: http://crd-legacy.lbl.gov/~xiaoye/SuperLU/#superlu_mt
  - requires KLU solver: http://faculty.cse.tamu.edu/davis/suitesparse.html
  - hint: add `-DSUPERLUMP_ENABLE=ON -DKLU_ENABLE=ON` to the `cmake` command to enable support

## set-up

General compilation parameters:
```
export CC=mpicc
export CXX=mpicxx
export WORKSPACE_PATH=/home/username/Workspace
```

HPX compilation parameters:
```
export HPX_INSTALL_PATH=$WORKSPACE_PATH/hpx-install
export PKG_CONFIG_PATH=$HPX_INSTALL_PATH/lib/pkgconfig/:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$HPX_INSTALL_PATH/lib:$LD_LIBRARY_PATH
export PATH=$HPX_INSTALL_PATH/bin:$PATH
```

neurox compilation parameters
```
export NEUROX_INSTALL_PATH=$WORKSPACE_PATH/neurox-install
export LD_LIBRARY_PATH=$NEUROX_INSTALL_PATH/lib:$LD_LIBRARY_PATH
export PATH=$NEUROX_INSTALL_PATH/bin:$PATH
```

## Compilation

```
cmake .. -DCMAKE_INSTALL_PREFIX=$NEUROX_INSTALL_PATH
```
*Note*: requires libcoreneuron installation with export of function pointers (`-DEXPORT_MECHS_FUNCTIONS`) and extra mechanisms path when applicable(`-DADDITIONAL_MECHPATH` and `-DADDITIONAL_MECHS`). To be on the safe side, we will also disable timing out of coreneuron after a certain period of inactivity (`-DDISABLE_NRN_TIMEOUT`) and OpenMP that may interfere with HPX scheduler (`-DCORENEURON_OPENMP`).
To allow header files to be exposed externally we activate the flag`ENABLE_DEV_FILES_INSTALLATION`. Vectorization is provided by setting the input data-structure as Array-Of-Structures or Structures-of-Arrays (`-DENABLE_SOA`). All other flags are optional and are added to speed-up compilation.
```
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_INSTALL_PREFIX=$NEUROX_INSTALL_PATH \
         -DADDITIONAL_MECHPATH=$NEURODAMUS_LIB_PATH \
         -DADDITIONAL_MECHS=$ADDITIONAL_MECHS_PATH \
         -DEXPORT_MECHS_FUNCTIONS=ON \
         -DCORENEURON_OPENMP=OFF \
         -DDISABLE_NRN_TIMEOUT=ON\
         -DENABLE_DEV_FILES_INSTALLATION=ON \
         -DCMAKE_CXX_FLAGS="-std=c++11" \  
         \
         -DENABLE_CUDA_MODULES=OFF \
         -DENABLE_NET_RECEIVE_BUFFERING=OFF \
         -DENABLE_OPENACC_INFO=OFF \
         -DENABLE_REPORTINGLIB=OFF \
         -DUNIT_TESTS=OFF \
         -DFUNCTIONAL_TESTS=OFF \
         -DCORENEURON_MAIN=OFF
```

## Execution

```
./neurox --help for execution parameters
./neurox -d <input-data-folder> -e <execution-time-milisecs> for execution with minimum parameters
```

## Misc

- We follow the google coding style (https://google.github.io/styleguide/cppguide.html) and format. To automatically format the code recursively in all folders use `clang-format`.
```
find ./neurox -iname *.h -o -iname *.cc  | xargs clang-format -i -style=Google
```

- To check for coding style errors, use `ccplint.py` (https://raw.githubusercontent.com/google/styleguide/gh-pages/cpplint/cpplint.py ). False positives can be ignored by putting `// NOLINT` at the end of the line or `// NOLINTNEXTLINE` in the previous line.
```
find ./neurox -iname *.h -o -iname *.cc  | xargs cpplint.py
```

- Documentation follows the doxygen (www.doxygen.org) notation and can be exported with `make doc`; 

- Developer note: to collect CPU hotspots with inter VTune (`amplxe-cl --help` for details):
```
amplxe-cl -collect hotspots -r r000hs srun -n 1 ./bin/neurox_exec -d /gpfs/bbp.cscs.ch/project/proj16/bmagalha/Circuits/PCP/1/coreneuron_input
```

## Copyright

- Blue Brain Project, École Polytechnique Fédérale de Lausanne (EPFL), Switzerland;
- Center for Research in Extreme Scale Technologies (CREST), Indiana University
