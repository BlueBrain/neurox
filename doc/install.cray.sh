set -x
module unload PrgEnv-pgi
module load PrgEnv-cray
module swap cce/8.4.0 cce/8.3.12 #on titan
module load craype-accel-nvidia35
module load perftools
module load cmake3/3.1.0
module load cudatoolkit

BASE=`pwd`
SRC_DIR=$BASE/sources
INSTALL_DIR=$BASE/install

#remember, -g -hprofile_generate disables function inling, why??
#COMPILATION_FLAG="-O3 -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP -hlist=a -h vector3 -hpragma=acc"
COMPILATION_FLAG=" -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP -hlist=a -h vector3 -hpragma=acc -hacc_model=auto_async_none"
COMPILATION_FLAG=" -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP -hlist=a -h vector3 -hpragma=acc -O3"
COMPILATION_FLAG=" -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP -hlist=a -h vector3 -hpragma=acc -hacc_model=auto_async_none -O2"

#for pgi -Mnodepchk might be problematic as it will disable all global checks
#module swap PrgEnv-cray PrgEnv-pgi
#COMPILATION_FLAG="-O3 -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP  -Minfo=all"

#directories where executables will be copied
GPU_EXE_DIR=$MEMBERWORK/csc192/gpu/
CPU_EXE_DIR=$MEMBERWORK/csc192/cpu/

#mechanism directory and list of files
MECHPATH=/ccs/home/kumbhar/workarena/systems/titan/repos/bbp/eurohack15/sources/neurodamus/lib/modlib
MECHS=/ccs/home/kumbhar/workarena/systems/titan/repos/bbp/eurohack15/sources/neurodamus/lib/modlib/coreneuron_modlist.txt

rm $GPU_EXE_DIR/coreneuron_exec
rm $CPU_EXE_DIR/coreneuron_exec

mkdir -p $SRC_DIR
mkdir -p $INSTALL_DIR

CORENEURON_SRC=$SRC_DIR/coreneuron
cd $CORENEURON_SRC

#mod2c path
export MODLUNIT=$SRC_DIR/mod2c/share/nrnunits.lib
export PATH=$PATH:$INSTALL_DIR/bin

export CC=`which cc`
export CXX=`which CC`

mkdir -p build_hpcopt_acc
cd build_hpcopt_acc

cmake .. -DMPI_C_LIBRARIES=$MPICH_DIR/lib  -DCORENEURON_LIBRARY_TYPE=STATIC  -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR -DMPI_C_INCLUDE_PATH=$MPICH_DIR/include -DMPI_C_COMPILER=`which cc` -DMPI_CXX_COMPILER=`which CC`  -DCMAKE_C_FLAGS="$COMPILATION_FLAG -D_CRAY_ACC_DEBUG"  -DCMAKE_CXX_FLAGS="$COMPILATION_FLAG" -DDISABLE_NRN_TIMEOUT=ON -DENABLE_SELECTIVE_GPU_PROFILING=ON -DADDITIONAL_MECHPATH=$MECHPATH -DADDITIONAL_MECHS=$MECHS -DCOMPILE_LIBRARY_TYPE=STATIC
#make clean
make VERBOSE=12 -j12
cd apps
rm ../bin/coreneuron_exec

cd $SRC_DIR/coreneuron/coreneuron/nrniv/
nvcc -c cuda_profile.cu
cd -
/opt/cray/craype/2.4.2/bin/CC   -O3 -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP -hlist=a -h vector3 -hpragma=acc -hacc_model=auto_async_none       CMakeFiles/coreneuron_exec.dir/main.cpp.o  -o ../bin/coreneuron_exec  ../coreneuron/libcoreneuron.a -Wc,$BASE/../cunrnran123/nrnran123.ptx $BASE/../cunrnran123/libnrnran123cuda.a $SRC_DIR/coreneuron/coreneuron/nrniv/cuda_profile.o

cd ..
echo `pwd`

mkdir -p $GPU_EXE_DIR
cp ./bin/coreneuron_exec $GPU_EXE_DIR/coreneuron_exec
cd $GPU_EXE_DIR

rm coreneuron_exec+* coreneuron_exec_*
pat_build -u -w coreneuron_exec -o coreneuron_exec_utrace
pat_build -O apa coreneuron_exec -o coreneuron_exec_osamp

#-u is for user but it skips all functions from acc as they are too small (?)

#CPU BUILD
COMPILATION_FLAG=" -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP -hlist=a -h vector3 -hpragma=acc -hnoacc -O2"
module unload craype-accel-nvidia35

cd $CORENEURON_SRC
mkdir -p build_hpcopt_x86
cd build_hpcopt_x86

cmake .. -DMPI_C_LIBRARIES=$MPICH_DIR/lib  -DCORENEURON_LIBRARY_TYPE=STATIC  -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR -DMPI_C_INCLUDE_PATH=$MPICH_DIR/include -DMPI_C_COMPILER=`which cc` -DMPI_CXX_COMPILER=`which CC`  -DCMAKE_C_FLAGS="$COMPILATION_FLAG"  -DCMAKE_CXX_FLAGS="$COMPILATION_FLAG" -DDISABLE_NRN_TIMEOUT=ON -DADDITIONAL_MECHPATH=$MECHPATH -DADDITIONAL_MECHS=$MECHS -DCOMPILE_LIBRARY_TYPE=STATIC
#make clean
rm ./bin/coreneuron_exec
make VERBOSE=12 -j12

mkdir -p $CPU_EXE_DIR
cp bin/coreneuron_exec $CPU_EXE_DIR/coreneuron_exec
cd $CPU_EXE_DIR
rm coreneuron_exec+* coreneuron_exec_*
pat_build -u coreneuron_exec -o coreneuron_exec_utrace
pat_build -O apa coreneuron_exec -o coreneuron_exec_osamp

#for DDT debugging
#/opt/cray/craype/2.4.2/bin/CC   -O3 -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP -hlist=a -h vector3 -hpragma=acc -hacc_model=auto_async_none       CMakeFiles/coreneuron_exec.dir/main.cpp.o  -o ../bin/coreneuron_exec  ../lib/libcoreneuron.a -Wc,/users/kumbhar/workarena/systems/daint/bbp/repos/curandom123/nrnran123.ptx /users/kumbhar/workarena/systems/daint/bbp/repos/curandom123/libnrnran123cuda.a -L/apps/daint/5.2.UP02/ddt/5.0.1-42607/lib/64 -Wl,--undefined=malloc -ldmallocthcxx

