# Makefile for Radex
# Uncomment the FC/FFLAGS combination of your choice
# Option 1:
#FC      = g77
#FFLAGS += -fno-automatic -Wall -O2
# Option 2:
FC =gfortran
f2pyFC=gnu95 #this must match FC (ifort=intelem, gfortran=gnu95)
FFLAGS += -g -fPIC
EXT_SUFFIX = $(shell python3-config --extension-suffix)

# Option 3:
#FC      = ifort
#FFLAGS += -O2
# Option 4:
#FC = g95
#FFLAGS += -Wall -O2
# there may be other options of course ...

OBJECTS = types.o slatec.o commondata.o solver.o io.o background.o  main.o
SOURCES = types.f90 slatec.f90 commondata.f90 solver.f90 io.f90 background.f90
EXEC    = radex
BINDIR  = ../bin/

all: $(EXEC) install clean

$(EXEC): $(OBJECTS) solver.o  Makefile main.f90
	$(LINK.F) $(OBJECTS) -o $@
	strip $@

python:  $(OBJECTS) solver.o  Makefile wrap.f90
	python3 -m numpy.f2py -c --f2cmap --fcompiler=${f2pyFC} $(OBJECTS) wrap.f90 -m radexwrap  -h wrap.pyf
	mv *${EXT_SUFFIX} ../radexwrap.so

%.o: %.f90 Makefile
	$(COMPILE.F) $*.f90 -o $@

install: $(EXEC)
	install -m 755 -p -s $(EXEC) $(BINDIR)

clean:
	- rm *.o *.mod *.so
	- rm $(EXEC)
