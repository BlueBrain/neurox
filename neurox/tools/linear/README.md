How to compile on BB5 with `likwid` flags:
```
$CXX -std=c++11 -DNDEBUG -DLIKWID_PERFMON -I/gpfs/bbp.cscs.ch/apps/tools/install/linux-rhel7-x86_64/gcc-4.8.5/likwid-4.3.0-vkur7z/include -L/gpfs/bbp.cscs.ch/apps/tools/install/linux-rhel7-x86_64/gcc-4.8.5/likwid-4.3.0-vkur7z/lib -llikwid ./benchmark.cc
```

How to run:
```
module load likwid
likwid-perfctr -C S0:0-0 -g CACHES -m  ./a.out 100 > nl0.txt
```
