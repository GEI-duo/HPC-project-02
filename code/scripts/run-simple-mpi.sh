#!/bin/bash

## Specifies the interpreting shell for the job.
#$ -S /bin/bash

## Execute the job from the current working directory.
#$ -cwd

MPICH_MACHINES=$TMPDIR/mpich_machines
cat $PE_HOSTFILE | awk '{print $1":"$2}' > $MPICH_MACHINES

export OMP_NUM_THREADS=$OMP_NUM_THREADS

## In this line you have to write the command that will execute your application.
mpiexec -f $MPICH_MACHINES -n $NSLOTS ./heat_mpi.o $SIZE $STEPS "tests/$MPI_PROGRAM"-"$MPI_NUM_JOBS"-"$STEPS"-"$SIZE".bmp


rm -rf $MPICH_MACHINES

