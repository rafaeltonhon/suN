FC=gfortran
FFLAGS=-c -O3
ldflags=-Wall -Wextra
sources=lattice.f90 mathtools.f90 thermalize.f90 measures.f90 suN.f90
objects=$(sources:.f90=.o)
executable=suN

all: $(sources) $(executable)
	$(FC) $(sources) $(FFLAGS) $(objects)

$(executable): $(objects)
	$(FC) $(ldflags) $(objects) -o $@

.f90.o:
	$(FC) $(FFLAGS) $< -o $@
