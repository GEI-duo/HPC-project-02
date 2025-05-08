#!/bin/bash

## Specifies the interpreting shell for the job.
#$ -S /bin/bash

## Specifies parallel environment
#$ -pe smp 4

## Execute the job from the current working directory.
#$ -cwd 

## In this line you have to write the command that will execute your application.
./"$OMP_PROGRAM".o $SIZE $STEPS "tests/$OMP_PROGRAM"-"$OMP_NUM_THREADS"-"$STEPS"-"$SIZE".bmp
