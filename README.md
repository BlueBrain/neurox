# neurox

A large-scale parallel+distributed asynchronous simulator of extended Hodgkin-Huxley neuron models

## Installation

### Dependencies
- cmake 2.8.12+: www.cmake.org
- hpx 4+: https://hpx.crest.iu.edu/
- tclap: http://tclap.sourceforge.net/
- mod2c: https://github.com/BlueBrain/mod2c
- libCoreNeuron: https://github.com/BlueBrain/coreneuron 

### set-up
```
#General compilation parameters
export CC=mpicc.mpich
export CXX=mpicxx.mpich
export WORKSPACE_PATH=/home/username/Workspace

#HPX compilation parameters
export HPX_INSTALL_PATH=$WORKSPACE_PATH/hpx-install
export PKG_CONFIG_PATH=$HPX_INSTALL_PATH/lib/pkgconfig/:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$HPX_INSTALL_PATH/lib:$LD_LIBRARY_PATH
export PATH=$HPX_INSTALL_PATH/bin:$PATH

#neurox compilation parameters
export NEUROX_INSTALL_PATH=$WORKSPACE_PATH/neurox-install
export MOD2C_INSTALL_PATH=$NEUROX_INSTALL_PATH
export LD_LIBRARY_PATH=$NEUROX_INSTALL_PATH/lib:$LD_LIBRARY_PATH
export PATH=$NEUROX_INSTALL_PATH/bin:$PATH
```

### Compilation with provided mechanisms
```
cmake .. -DCMAKE_INSTALL_PREFIX=$NEUROX_INSTALL_PATH  #add-DENABLE_MPI=OFF to buid with SMP
```

### Compilation with user-specified mechanisms
```
export NEURODAMUS_LIB_PATH=$WORKSPACE_PATH/neurodamus/lib/modlib
export ADDITIONAL_MECHS_PATH=$WORKSPACE_PATH/neurodamus/lib/modlib/coreneuron_modlist.txt

cmake .. -DADDITIONAL_MECHPATH=$NEURODAMUS_LIB_PATH \
         -DADDITIONAL_MECHS=$ADDITIONAL_MECHS_PATH  \
         -DCMAKE_INSTALL_PREFIX=$NEUROX_INSTALL_PATH \
         #add -DENABLE_MPI=OFF to build with SMP
          
```

## Execution
- `./neurox --help` for execution parameters
- `./neurox -d <input-data-folder> -e <execution-time-milisecs>` for execution with minimum parameters

## Copyright 
Blue Brain Project, EPFL, Switzerland; and Center for Research in Extreme Scale Technologies, Indiana University
