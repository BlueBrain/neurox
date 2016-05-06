#! /bin/sh 

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/bmagalha/Workspace/neurox/build/lib
export OMP_NUM_THREADS=1

# Run the executable
#SRUN_EXTRA="--propagate=CORE"
SRUN_EXTRA=
if [ -n "$VALGRIND" -a -n "$VALGRIND_PRELOAD" ]; then
    echo "Running with valgrind"
    LD_PRELOAD=$VALGRIND_PRELOAD \
    /usr/bin/srun -n 1 $SRUN_EXTRA $VALGRIND /home/bmagalha/Workspace/neurox/build/bin/coreneuron_exec --datpath=/home/bmagalha/Workspace/neurox/tests/integration/ring_IClamp --outpath=/home/bmagalha/Workspace/neurox/build/tests/integration/ring_IClamp --celsius=6.3 --tstop=100. -mpi 
else
    /usr/bin/srun -n 1 $SRUN_EXTRA /home/bmagalha/Workspace/neurox/build/bin/coreneuron_exec --datpath=/home/bmagalha/Workspace/neurox/tests/integration/ring_IClamp --outpath=/home/bmagalha/Workspace/neurox/build/tests/integration/ring_IClamp --celsius=6.3 --tstop=100. -mpi 
fi
exitvalue=$?

# Check for error result
if [ $exitvalue -ne 0 ]; then
  echo "Error status value: $exitvalue"
  exit $exitvalue
fi

# diff outputed files with reference
cd /home/bmagalha/Workspace/neurox/build/tests/integration/ring_IClamp

if [ ! -f out0.dat ]
then
  echo "No output files. Test failed!"
  exit 1
fi

cat out[0-9]*.dat > out_cn.dat
sort -k 1n,1n -k 2n,2n out_cn.dat > sort_out.dat
diff sort_out.dat /home/bmagalha/Workspace/neurox/tests/integration/ring_IClamp/out.dat.ref > diff.dat

if [ -s diff.dat ] 
then
  echo "Results are different, check the file diff.dat. Test failed!"
  exit 1
else
  echo "Results are the same, test passed"
  rm -f *.dat
  exit 0
fi
