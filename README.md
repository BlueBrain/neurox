# neurox

A large-scale parallel+distributed asynchronous simulator of extended Hodgkin-Huxley neuron models

## Installation

### Dependencies
- cmake 2.8.12+: www.cmake.org
- hpx 4+: https://hpx.crest.iu.edu/
- tclap: http://tclap.sourceforge.net/
- libCoreNeuron: https://github.com/brunomaga/coreneuron 

### set-up
```
#General compilation parameters
export CC=mpicc
export CXX=mpicxx
export WORKSPACE_PATH=/home/username/Workspace

#HPX compilation parameters
export HPX_INSTALL_PATH=$WORKSPACE_PATH/hpx-install
export PKG_CONFIG_PATH=$HPX_INSTALL_PATH/lib/pkgconfig/:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$HPX_INSTALL_PATH/lib:$LD_LIBRARY_PATH
export PATH=$HPX_INSTALL_PATH/bin:$PATH

#neurox compilation parameters
export NEUROX_INSTALL_PATH=$WORKSPACE_PATH/neurox-install
export LD_LIBRARY_PATH=$NEUROX_INSTALL_PATH/lib:$LD_LIBRARY_PATH
export PATH=$NEUROX_INSTALL_PATH/bin:$PATH
```

### Compilation
```
cmake .. -DCMAKE_INSTALL_PREFIX=$NEUROX_INSTALL_PATH
```
NOTE, requires libcoreneuron installation with export of function pointers, and mechanisms functions e.g.:
```
cmake .. -DCMAKE_INSTALL_PREFIX=$NEUROX_INSTALL_PATH \
         -DADDITIONAL_MECHPATH=$NEURODAMUS_LIB_PATH \
         -DADDITIONAL_MECHS=$ADDITIONAL_MECHS_PATH \
         -DEXPORT_MECHS_FUNCTIONS=ON
```

## Execution
- `./neurox --help` for execution parameters
- `./neurox -d <input-data-folder> -e <execution-time-milisecs>` for execution with minimum parameters

## Copyright 
Blue Brain Project, EPFL, Switzerland; and Center for Research in Extreme Scale Technologies, Indiana University
