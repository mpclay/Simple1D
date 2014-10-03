FC=gfortran
FFLAGS= -O3

simple1D: Simple1D.F90 parameters.o boundaryconditions.o flux.o
	${FC} ${FFLAGS} Simple1D.F90 Parameters.o BoundaryConditions.o Flux.o -o simple1D.x

flux.o: parameters.o
	${FC} ${FFLAGS} -c Flux.F90

boundaryconditions.o: parameters.o
	${FC} ${FFLAGS} -c BoundaryConditions.F90

parameters.o:
	${FC} ${FFLAGS} -c Parameters.F90

clean:
	rm -rf *.o *.mod *.dat

veryclean: clean
	rm -rf simple1D.x

