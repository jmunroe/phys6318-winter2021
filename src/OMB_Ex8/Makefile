FC=gfortran
FFLAGS=-g -O3 

run.exe: main.o sub.o param.o io.o
	$(FC) $(FFLAGS) -o run.exe main.o sub.o param.o io.o -L$(CONDA_PREFIX)/lib -lnetcdff

main.o: main.f95 param.o sub.o io.o
	$(FC) $(FFLAGS) -g -c main.f95

sub.o: sub.f95
	$(FC) $(FFLAGS) -g -c sub.f95

param.o: param.f95
	$(FC) $(FFLAGS) -g -c param.f95

io.o: io.f95 param.o
	$(FC) $(FFLAGS) -g -c io.f95 -I $(CONDA_PREFIX)/include

.PHONY: clean
clean:
	rm *.o
	rm *.mod
