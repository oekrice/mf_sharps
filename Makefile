# Set compiler, flags and netcdf library according to machine we are using:



TARGET = mf3d
$(info Machine name: $(shell hostname))

MPIF90 ?= /usr/lib64/openmpi/bin/mpifort
NETCDF = -I /usr/lib64/gfortran/modules
NETCDFLIB = -L/usr/lib64/libnetcdff.so.7 -lnetcdff
FFLAGS = -O3 #-fcheck=all -Wuninitialized -march=native -fimplicit-none -Wall -Wextra -ffast-math -funroll-loops --param max-unroll-times=5

# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------
SRCFILES = shared_data.o mpi_tools.o pressure.o boundary.o evolve.o init.o output.o main.o

$(info Compiling (or uncompiling)...)

all: main

SDF := SDF/FORTRAN
SDFMOD = $(SDF)/include/sdf.mod
SRCDIR = src
OBJDIR = obj
BINDIR = bin
FC = $(MPIF90)
DATE := $(shell date +%s)
MACHINE := $(shell uname -n)
  $(D)_MACHINE='"$(MACHINE)"'


OBJFILES := $(SRCFILES:.f90=.o)
OBJFILES := $(OBJFILES:.F90=.o)

FULLTARGET = $(BINDIR)/$(TARGET)

VPATH = $(SRCDIR):$(SRCDIR)/core:$(SDF)/src:$(OBJDIR)

-include $(SRCDIR)/COMMIT

ifeq ($(DONE_COMMIT),)
main: commit
else
main: $(FULLTARGET)
endif

# Rule to build the fortran files

%.o: %.f90
	@mkdir -p $(BINDIR) $(OBJDIR)
	$(FC) -c $(FFLAGS) $(NETCDF) -J $(OBJDIR) -o $(OBJDIR)/$@ $<

%.o: %.F90
	@mkdir -p $(BINDIR) $(OBJDIR)
	$(FC) -c $(FFLAGS) $(NETCDF) -J $(OBJDIR) -o $(OBJDIR)/$@ $<

$(FULLTARGET): $(OBJFILES)
	$(FC) $(FFLAGS) -J $(OBJDIR) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES)) $(NETCDFLIB)


clean:
	@rm -rf $(BINDIR) $(OBJDIR)


#$(OBJFILES): | $(OBJDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

commit: FORCE
	@sh $(SRCDIR)/gen_commit_string.sh && $(MAKE) $(MAKECMDGOALS) DONE_COMMIT=1

FORCE:

.PHONY: commit clean cleanall tidy datatidy visit visitclean main FORCE

shared_data.o: shared_data.f90
mpi_tools.o: mpi_tools.f90
pressure.o: pressure.f90 shared_data.o
evolve.o: evolve.f90 shared_data.o mpi_tools.o
init.o: init.f90 shared_data.o mpi_tools.o
output.o: output.f90 shared_data.o
main.o: main.f90 shared_data.o init.o evolve.o
boundary.o: boundary.f90 shared_data.o mpi_tools.o
