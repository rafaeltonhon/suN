FC=gfortran
FFLAGS=-O3 -g -fcheck=all

.SUFFIXES:
.SUFFIXES: .o .f90

.f90.o:
	$(FC) $(FFLAGS) -c $<

OBJECTS=\
   suN.o\
   lattice.o\
   mathtools.o\
   thermalize.o\
   measures.o\

SRC=\
   suN.f90\
   lattice.f90\
   mathtools.f90\
   thermalize.f90\
   measures.f90\

suN: $(OBJECTS)
	$(FC) -o $@ $^

clean:
	rm *\.o suN

ctags: $(SRC)
	ctags $(SRC)

suN.o: lattice.o mathtools.o thermalize.o measures.o suN.f90
	$(FC) -Wall -c suN.f90 -o suN.o

thermalize.o: lattice.o mathtools.o thermalize.f90
	$(FC) -Wall -c thermalize.f90 -o thermalize.o

measures.o: lattice.o mathtools.o measures.f90
	$(FC) -Wall -c measures.f90 -o measures.o

mathtools.o: lattice.o mathtools.f90
	$(FC) -Wall -c mathtools.f90 -o mathtools.o
	
lattice.o: lattice.f90
	$(FC) -Wall -c lattice.f90 -o lattice.o

