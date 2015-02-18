EXECS=totalViewshed
MPICC?=/Users/atodesco/Documents/gis_algorithms/openMPI/bin/mpic++
#MPICC?= <your path to open mpi>/bin/mpic++

all: ${EXECS}

totalViewshed: totalViewshed.c  totalViewshed.h
	${MPICC} -o totalViewshed totalViewshed.c 

clean: 
	rm ${EXECS}
