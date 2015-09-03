set -x
#set -e

module load PrgEnv-cray
module load craype-accel-nvidia35
module load perftools
module load cmake

BASE=/scratch/daint/kumbhar/workspace
BASE=`pwd`
SRC_DIR=$BASE/sources
INSTALL_DIR=$BASE/install


COMPILATION_FLAG="-O3 -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP -hlist=a -h vector3"
#remember, -g -hprofile_generate disables function inling, why??
COMPILATION_FLAG="-O3 -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP -hlist=a -h vector3 -hpragma=acc -hacc_model=auto_async_none"

#for pgi -Mnodepchk might be problematic as it will disable all global checks
#module swap PrgEnv-cray PrgEnv-pgi
#COMPILATION_FLAG="-O3 -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP  -Minfo=all"

mkdir -p $SRC_DIR
mkdir -p $INSTALL_DIR

CORENEURON_SRC=$SRC_DIR/coreneuron
cd $CORENEURON_SRC

export MODLUNIT=$SRC_DIR/mod2c/share/nrnunits.lib
export PATH=$PATH:$INSTALL_DIR/bin

sed -i 's/-Wno-error//' $CORENEURON_SRC/coreneuron/CMakeLists.txt
sed -i 's/corebluron/sandbox\/kumbhar\/hpcopt_stochkv_tmp/' $CORENEURON_SRC/.gitexternals
sed -i 's/-std=c++11//' CMake/common/Compiler.cmake

export CC=`which cc`
export CXX=`which CC`

mkdir -p build_hpcopt_acc
cd build_hpcopt_acc

cmake .. -DMPI_C_LIBRARIES=$MPICH_DIR/lib  -DCORENEURON_LIBRARY_TYPE=STATIC  -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR -DMPI_C_INCLUDE_PATH=$MPICH_DIR/include -DMPI_C_COMPILER=`which cc` -DMPI_CXX_COMPILER=`which CC`  -DCMAKE_C_FLAGS="$COMPILATION_FLAG"  -DCMAKE_CXX_FLAGS="$COMPILATION_FLAG" -DDISABLE_NRN_TIMEOUT=ON
make VERBOSE=1 -j1

#for DDT debugging
#/opt/cray/craype/2.4.0/bin/CC   -O3 -g -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP -hlist=a -h vector3 -hpragma=acc -hacc_model=auto_async_none -hnoacc       CMakeFiles/coreneuron_exec.dir/main.cpp.o  -o ../bin/coreneuron_exec  ../lib/libcoreneuron.a -L/apps/daint/5.2.UP02/ddt/5.0.1-42607/lib/64 -Wl,--undefined=malloc -ldmallocthcxx

#to link with CUDA
cd /users/kumbhar/workarena/systems/daint/bbp/repos/hackathon15/sources/coreneuron/build_hpcopt_acc/apps
/opt/cray/craype/2.4.0/bin/CC   -O3 -DSWAP_ENDIAN_DISABLE_ASM -DLAYOUT=0 -DDISABLE_HOC_EXP -hlist=a -h vector3 -hpragma=acc -hacc_model=auto_async_none       CMakeFiles/coreneuron_exec.dir/main.cpp.o  -o ../bin/coreneuron_exec  ../lib/libcoreneuron.a -Wc,/users/kumbhar/workarena/systems/daint/bbp/repos/curandom123/nrnran123.ptx /users/kumbhar/workarena/systems/daint/bbp/repos/curandom123/libnrnran123cuda.a
