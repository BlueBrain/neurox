How to compile on BB5 with `likwid` flags:
```
$CXX -std=c++11 -DNDEBUG -DLIKWID_PERFMON -I/gpfs/bbp.cscs.ch/apps/tools/install/linux-rhel7-x86_64/gcc-4.8.5/likwid-4.3.0-vkur7z/include -L/gpfs/bbp.cscs.ch/apps/tools/install/linux-rhel7-x86_64/gcc-4.8.5/likwid-4.3.0-vkur7z/lib -llikwid ./benchmark.cc
```

How to run:
```
module load likwid
likwid-perfctr -C S0:0-0 -g CACHES -m  ./a.out 100 > nl0.txt
```


Bath script for slurm job submission:
```
#!/bin/bash
#SBATCH --job-name="cache-efficiency"
#SBATCH --time=8:00:00
#SBATCH --partition=prod
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --account=proj16
#SBATCH --mail-user=bruno.magalhaes@epfl.ch
#SBATCH --mail-type=ALL
#SBATCH --output=cache-eff.log
#SBATCH --error=cache-eff.err

NODES=1
L1_CACHE_SIZE=512

module load likwid

for SIZE in 2048 #16 32 64 128 256 512 1024 2048 4096
do
  echo "============= RUNNING SIZE=$SIZE ================="
  srun -n 1 -c 1 likwid-perfctr -C S0:0-0 -g CACHES -m ./nl0 $SIZE > cache_${SIZE}_nl0.txt
  srun -n 1 -c 1 likwid-perfctr -C S0:0-0 -g CACHES -m ./lin0 $SIZE > cache_${SIZE}_lin0.txt
#  srun -n 1 -c 1 likwid-perfctr -C S0:0-0 -g CACHES -m ./nl1 $SIZE > cache_${SIZE}_nl1.txt
#  srun -n 1 -c 1 likwid-perfctr -C S0:0-0 -g CACHES -m ./lin1 $SIZE > cache_${SIZE}_lin1.txt
done


for SIZE in 2048 #in 16 32 64 128 256 512 1024 2048 4096
do
  echo "============= RUNNING SIZE=$SIZE ================="
  srun -n 1 -c 1 likwid-perfctr -C S0:0-0 -g CACHES -m ./nl0s $SIZE > cache_${SIZE}_nl0s.txt
  srun -n 1 -c 1 likwid-perfctr -C S0:0-0 -g CACHES -m ./lin0s $SIZE > cache_${SIZE}_lin0s.txt
  #srun -n 1 -c 1 likwid-perfctr -C S0:0-0 -g CACHES -m ./nl1s $SIZE > cache_${SIZE}_nl1s.txt
  #srun -n 1 -c 1 likwid-perfctr -C S0:0-0 -g CACHES -m ./lin1s $SIZE > cache_${SIZE}_lin1s.txt
done
```

where `nl0` is the std based data strctures benchmark (disabled flag `#define LINEAR` in `benchmark.cc`), `lin0` is the linear based data structures, and binaries suffixed with `*s` enable flag `SCHEDULER`.
