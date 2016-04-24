module unload PrgEnv-pgi
module load PrgEnv-cray
module swap cce/8.4.0 cce/8.3.12 #on titan
module load craype-accel-nvidia35
module load cce/8.3.12
module load perftools

#module load cmake/2.8.11.2
#module unload cray-libsci_acc/3.3.0
#module unload cudatoolkit
#module load cudatoolkit/6.5.14-1.0502.9836.8.1
#module swap cudatoolkit/7.0.28-1.0502.10280.4.1 cudatoolkit/6.5.14-1.0502.9836.8.1
#module load cudatoolkit
#module swap cudatoolkit/7.0.28-1.0502.10280.4.1 cudatoolkit/6.5.14-1.0502.9836.8.1

#CMAKE issue where -no-shared or -no-fic added?
export PATH=/ccs/home/kumbhar/workarena/systems/titan/softwares/install/common/cmake-2.8.12.2/bin/:$PATH

#installation paths
BASE=`pwd`
SRC_DIR=$BASE/sources
INSTALL_DIR=$BASE/install

#Cray asynchronous launching: kernel, all or none
ACC_ASYNC_MODEL="auto_async_none"
ACC_ASYNC_MODEL="auto_async_all"

#Cray compielr options for floating optimizations
FP_OPTS=""
FP_OPTS="-h fp0"

#remember, -g -hprofile_generate disables function inling, why??
COMPILATION_FLAG=" -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP -hlist=a -h vector3 -hpragma=acc -O3"
COMPILATION_FLAG=" -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP -hlist=a -h vector3 -hpragma=acc -hacc_model=${ACC_ASYNC_MODEL} -O2 ${FP_OPTS}"

#for pgi -Mnodepchk might be problematic as it will disable all global checks
#module swap PrgEnv-cray PrgEnv-pgi
#COMPILATION_FLAG="-O3 -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP  -Minfo=all"

#directories where executables will be copied
GPU_EXE_DIR=$MEMBERWORK/csc192/gpu/
CPU_EXE_DIR=$MEMBERWORK/csc192/cpu/

#mechanism directory and list of files
MECHPATH=/ccs/home/kumbhar/workarena/systems/titan/repos/bbp/eurohack15/sources/neurodamus/lib/modlib
MECHS=/ccs/home/kumbhar/workarena/systems/titan/repos/bbp/eurohack15/sources/neurodamus/lib/modlib/coreneuron_modlist.txt

#remove executables
rm $GPU_EXE_DIR/coreneuron_exec
rm $CPU_EXE_DIR/coreneuron_exec

#create directories
mkdir -p $SRC_DIR
mkdir -p $INSTALL_DIR

CORENEURON_SRC=$SRC_DIR/coreneuron
cd $CORENEURON_SRC

#mod2c path
export MODLUNIT=$SRC_DIR/mod2c/share/nrnunits.lib
export PATH=$PATH:$INSTALL_DIR/bin

export CC=`which cc`
export CXX=`which CC`

#OpenACC build directory
mkdir -p build_hpcopt_acc
cd build_hpcopt_acc

#build OpenACC version
cmake .. -DMPI_C_LIBRARIES=$MPICH_DIR/lib  -DCORENEURON_LIBRARY_TYPE=STATIC  -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR -DMPI_C_INCLUDE_PATH=$MPICH_DIR/include -DMPI_CXX_INCLUDE_PATH=$MPICH_DIR/include -DMPI_C_COMPILER=`which cc` -DMPI_CXX_COMPILER=`which CC`  -DCMAKE_C_FLAGS="$COMPILATION_FLAG -D_CRAY_ACC_DEBUG"  -DCMAKE_CXX_FLAGS="$COMPILATION_FLAG" -DDISABLE_NRN_TIMEOUT=ON -DENABLE_SELECTIVE_GPU_PROFILING=ON -DADDITIONAL_MECHPATH=$MECHPATH -DADDITIONAL_MECHS=$MECHS -DCOMPILE_LIBRARY_TYPE=STATIC

#make clean
make VERBOSE=1 -j12

#todo: CMake changes for for .cu code
cd apps
rm ../bin/coreneuron_exec
cd $SRC_DIR/coreneuron/coreneuron/nrniv/
nvcc -c cuda_profile.cu

#link with random123 lib
cd -
/opt/cray/craype/2.4.2/bin/CC   -O3 -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP -hlist=a -h vector3 -hpragma=acc -hacc_model=${ACC_ASYNC_MODEL}       CMakeFiles/coreneuron_exec.dir/main.cpp.o  -o ../bin/coreneuron_exec  ../coreneuron/libcoreneuron.a -Wc,$BASE/../cunrnran123/nrnran123.ptx $BASE/../cunrnran123/libnrnran123cuda.a $SRC_DIR/coreneuron/coreneuron/nrniv/cuda_profile.o

cd ..

#copy executable
mkdir -p $GPU_EXE_DIR
cp ./bin/coreneuron_exec $GPU_EXE_DIR/coreneuron_exec
cd $GPU_EXE_DIR

#instrument with CrayPat, -u is for user defined functions but it skips all functions from acc as they might be too small (?)
#rm coreneuron_exec+* coreneuron_exec_*
pat_build -u -w coreneuron_exec -f -o coreneuron_exec_utrace
pat_build -O apa coreneuron_exec -f -o coreneuron_exec_osamp


#---------------------------CUDA VERSION COMPILATION---------------------------------
if true;
then
    #cuda version
    cd $CORENEURON_SRC
    mkdir -p build_hpcopt_acc_cuda
    cd build_hpcopt_acc_cuda

    cmake .. -DMPI_C_LIBRARIES=$MPICH_DIR/lib  -DCORENEURON_LIBRARY_TYPE=STATIC  -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR -DMPI_C_INCLUDE_PATH=$MPICH_DIR/include -DMPI_C_COMPILER=`which cc` -DMPI_CXX_COMPILER=`which CC`  -DCMAKE_C_FLAGS="$COMPILATION_FLAG -D_CRAY_ACC_DEBUG -DENABLE_CUDA -DENABLE_CUDA_INTERFACE"  -DCMAKE_CXX_FLAGS="$COMPILATION_FLAG -DENABLE_CUDA -DENABLE_CUDA_INTERFACE" -DDISABLE_NRN_TIMEOUT=ON -DENABLE_SELECTIVE_GPU_PROFILING=ON -DADDITIONAL_MECHPATH=$MECHPATH -DADDITIONAL_MECHS=$MECHS -DCOMPILE_LIBRARY_TYPE=STATIC

    #make clean
    make VERBOSE=1 -j12

    cd apps

    #rm ../bin/coreneuron_exec
    cd ~/workarena/systems/titan/repos/bbp/eurohack15/cuda

    #build add cuda kernels
    make -j12 && cd -

    #link against cuda kernels and build library
    /opt/cray/craype/2.4.2/bin/CC   -O3 -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP -hlist=a -h vector3 -hpragma=acc -hacc_model=${ACC_ASYNC_MODEL}       CMakeFiles/coreneuron_exec.dir/main.cpp.o  -o ../bin/coreneuron_exec  ../coreneuron/libcoreneuron.a -Wc,$BASE/../cunrnran123/nrnran123.ptx $SRC_DIR/coreneuron/coreneuron/nrniv/cuda_profile.o /ccs/home/kumbhar/workarena/systems/titan/repos/bbp/eurohack15/cuda/*.o
    cp ../bin/coreneuron_exec $GPU_EXE_DIR/cu_coreneuron_exec
fi

########----------------------CPU Build------------------------------------

#CPU BUILD
#COMPILATION_FLAG=" -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP -O2 -opt-report5"
COMPILATION_FLAG=" -DENABLE_SELECTIVE_PROFILING -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP -hlist=a -h vector3 -hpragma=acc -hnoacc -O2 ${FP_OPTS}"
module unload craype-accel-nvidia35
#module unload perftools
#module load tau
#module swap PrgEnv-cray PrgEnv-intel
#module load gcc
export PATH=/ccs/home/kumbhar/workarena/systems/titan/softwares/:$PATH
#export TAU_MAKEFILE=/sw/xk6/tau/2.24.1/cle5.2_cray8.4.0/craycnl/lib/Makefile.tau-cray-papi-mpi-pdt-openmp-opari
#export TAU_OPTIONS="-optRevert -optTauSelectFile=/ccs/home/kumbhar/tau.instrument"
#export CC=cc.tau
#export CXX=CC.tau

cd $CORENEURON_SRC
mkdir -p build_hpcopt_x86
cd build_hpcopt_x86

cmake .. -DMPI_C_LIBRARIES=$MPICH_DIR/lib  -DCORENEURON_LIBRARY_TYPE=STATIC  -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR -DMPI_C_INCLUDE_PATH=$MPICH_DIR/include -DMPI_C_COMPILER=`which cc` -DMPI_CXX_COMPILER=`which CC`  -DCMAKE_C_FLAGS="$COMPILATION_FLAG"  -DCMAKE_CXX_FLAGS="$COMPILATION_FLAG" -DDISABLE_NRN_TIMEOUT=ON -DADDITIONAL_MECHPATH=$MECHPATH -DADDITIONAL_MECHS=$MECHS -DCOMPILE_LIBRARY_TYPE=STATIC
#make clean
#rm -f ./bin/coreneuron_exec
make VERBOSE=1 -j12

#copy executable
mkdir -p $CPU_EXE_DIR
cp bin/coreneuron_exec $CPU_EXE_DIR/coreneuron_exec

#instrument with CrayPat
cd $CPU_EXE_DIR
#rm coreneuron_exec+* coreneuron_exec_*
pat_build -u coreneuron_exec -f -o coreneuron_exec_utrace
pat_build -O apa coreneuron_exec -f -o coreneuron_exec_osamp

#for DDT debugging
#/opt/cray/craype/2.4.2/bin/CC   -O3 -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP -hlist=a -h vector3 -hpragma=acc -hacc_model=${ACC_ASYNC_MODEL}       CMakeFiles/coreneuron_exec.dir/main.cpp.o  -o ../bin/coreneuron_exec  ../lib/libcoreneuron.a -Wc,/users/kumbhar/workarena/systems/daint/bbp/repos/curandom123/nrnran123.ptx /users/kumbhar/workarena/systems/daint/bbp/repos/curandom123/libnrnran123cuda.a -L/apps/daint/5.2.UP02/ddt/5.0.1-42607/lib/64 -Wl,--undefined=malloc -ldmallocthcxx

