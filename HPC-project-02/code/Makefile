# Compiler and flags
CC = mpicc
CFLAGS = -fopenmp
SCRIPTS=./scripts

# Source files
SOURCES = $(wildcard *.c)

# Object files
OBJECTS = $(SOURCES:.c=.o)

# Rule to compile object files
%.o: %.c
	$(CC) $(CFLAGS) $< -o $@

# Compile only target
compile: $(OBJECTS)


serial: compile
	${SCRIPTS}/run-multiple.sh serial

mpi: compile
	${SCRIPTS}/run-multiple.sh mpi

all: compile mpi

# Extract test time results
time:
	${SCRIPTS}/run-extract-results.sh

# Clean target
clogs:
	rm -f omp_*.e* omp_*.o* mpi_*.e* mpi_*.o* serial_*.e* serial_*.o*

ctests:
	rm -f tests/* 

ctime:
	rm -f results.csv

clean: ctests clogs ctime
	rm -f $(OBJECTS)

test:
	${SCRIPTS}/run-tests.sh

# Phony targets to prevent conflicts with files named 'compile' or 'clean'
.PHONY: compile clean
