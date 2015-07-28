rm out*.*
export OMP_NUM_THREADS=1

module load PrgEnv-cray
module load craype-accel-nvidia35
module load perftools
module load cmake

aprun -N 1 -n 1 --cc none  -d 2 /scratch/daint/kumbhar/workspace/exe/coreneuron-soa-daint.cray.cn-master.neuroda-sandbox_kumbhar_hpcopt.prof-none-vec -mpi -e 15 -d /scratch/daint/kumbhar/data/hpcopt/4k/64/ -f files.dat
sort -k 1n,1n -k 2n,2n out0.dat > out0.dat.master.sorted

rm out0.dat 

aprun -N 1 -n 1 --cc none  -d 2 /users/kumbhar/workarena/systems/daint/bbp/repos/hackathon15/sources/coreneuron-old/build_hpcopt_x86/bin/coreneuron_exec -mpi -e 15 -d /scratch/daint/kumbhar/data/hpcopt/4k/64/ -f files.dat
sort -k 1n,1n -k 2n,2n out0.dat > out0.dat.acc-base.sorted

rm out0.dat 


aprun -N 1 -n 1 --cc none  -d 2 -b nvprof -o gpu.nvprof /users/kumbhar/workarena/systems/daint/bbp/repos/hackathon15/sources/coreneuron/build_hpcopt_acc/bin/coreneuron_exec -mpi -e 15 -d /scratch/daint/kumbhar/data/hpcopt/4k/64/ -f files.dat

sleep 4

sort -k 1n,1n -k 2n,2n out0.dat > out0.dat.gpu.sorted
