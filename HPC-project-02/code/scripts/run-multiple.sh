#!/bin/bash

prog=$1

case $prog in
    "serial" )
        program="heat_serial"
        ;;
    "omp" )
        program="heat_omp"
        ;;
    "mpi" )
	program="heat_mpi"
	;;
     * )
        echo "Unknown program: $program. Options are 'serial', 'omp' and 'mpi'."
        exit 1
esac

threads=(1 2)
threads=(1 2 4 8 16 32)
steps=(100 1000)
steps=(100 1000 10000 100000)
matrix=(100 1000)
matrix=(100 1000 2000)

for t in "${threads[@]}"; do
    for s in "${steps[@]}"; do
        for m in "${matrix[@]}"; do
		if [ "${program}" = "heat_mpi" ]; then
			echo -e "===Running program '$program' with {jobs=$t, steps=$s, matrix_size=$m}..."
			qsub -V -pe mpich "$t" -N "$prog"_"$t"_"$s"_"$m" -v MPI_NUM_JOBS="$t",OMP_NUM_THREADS="4",MPI_PROGRAM="$program",STEPS="$s",SIZE="$m" ./scripts/run-simple-mpi.sh	
			continue		
fi
            echo -e "===Running program '$program' with {threads=$t, steps=$s, matrix_size=$m}..."
            qsub -V -N "$prog"_"$t"_"$s"_"$m" -v OMP_NUM_THREADS="$t",OMP_PROGRAM="$program",STEPS="$s",SIZE="$m" ./scripts/run-simple-omp.sh
        done
    done
done

