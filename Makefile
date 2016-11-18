
# The header file provided:
#HEADERS = heat_serial.h 

# The compiler: gcc for C program:
CC = gcc

# Extra compiler flags: -Wall turns on most compiler warnings; -I sets path (here):
CFLAGS = -Wall -I./

# Makefile recipes:
all : heat_serial heat_omp heat_mpi

# Building the main drivers with integrator dependencies:
heat_serial: heat_serial.c heat_serial.o 
	$(CC) $(CFLAGS) -o heat_serial heat_serial.c -lm

heat_omp: heat_omp.c heat_omp.o 
	$(CC) $(CFLAGS) -fopenmp -o heat_omp heat_omp.c -lm

heat_mpi: heat_mpi.c heat_mpi.o 
	mpicc -o heat_mpi heat_mpi.c -lm
	#mpicc -c heat_mpi.c -o heat_mpi.o


# Implementations of the integrators:
heat_serial.o: heat_serial.c 
	$(CC) $(CFLAGS) -c heat_serial.c

heat_omp.o: heat_omp.c 
	$(CC) $(CFLAGS) -fopenmp -c heat_omp.c

heat_mpi.o: heat_mpi.c 
	mpicc -c heat_mpi.c


# Clean command to delete all relevant files:
clean: 
	rm heat_serial.o heat_omp.o heat_mpi.o
	rm heat_serial heat_omp heat_mpi
