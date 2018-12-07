How to compile (BB5 with likwid)
module load likwid
$CXX -std=c++11 -DLIKWID_PERFMON -I/gpfs/bbp.cscs.ch/apps/tools/install/linux-rhel7-x86_64/gcc-4.8.5/likwid-4.3.0-vkur7z/include -L/gpfs/bbp.cscs.ch/apps/tools/install/linux-rhel7-x86_64/gcc-4.8.5/likwid-4.3.0-vkur7z/lib -llikwid ./benchmark.cc

How to run:
likwid-perfctr -C S0:0-0 -g CACHES -m  ./a.out 100 > nl0.txt
