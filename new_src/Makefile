# GNU Fortran / MPI build for the current free-form radial solver.
#
# Default:  make -j10
# Checked:  make -j10 debug
# Custom:   make -j10 EXTRA_FFLAGS="..." EXTRA_LDFLAGS="..."

MPIFC ?= mpifort
EXEC  ?= qcv

# -fdefault-real-8 is still required by the not-yet-kind-explicit new_src
# fields and must agree with MPIFLOAT=MPI_DOUBLE_PRECISION.  The companion
# double flag prevents DOUBLE PRECISION from being promoted to 16 bytes.
COMMON_FFLAGS = -std=f2018 -fdefault-real-8 -fdefault-double-8 \
	-fimplicit-none -Wall -Wextra -Wimplicit-interface

# Keep the production build numerically conservative: no fast-math and no
# host-specific instruction set.  Architecture/offload work will be enabled
# explicitly after its numerical regression tests exist.
RELEASE_FFLAGS = -O2

# Runtime diagnostics for development.  Bounds checking is intentionally not
# present in the production flags because it is expensive in the inner loops.
DEBUG_FFLAGS = -O0 -g -fcheck=all -fbacktrace \
	-ffpe-trap=invalid,zero,overflow -finit-real=snan

FFLAGS = $(COMMON_FFLAGS) $(RELEASE_FFLAGS) $(EXTRA_FFLAGS)
LDFLAGS += $(EXTRA_LDFLAGS)
LDLIBS  += $(EXTRA_LDLIBS)

OBJ = he3_kinds.o legacy_anderson_mixing.o global_dec.o bulkgap.o \
	init_calc.o mpicalls.o interpol.o riccati.o getnewses.o packses.o \
	iter_AA.o iter.o qcv.o

.DELETE_ON_ERROR:

all: $(EXEC)

$(EXEC): $(OBJ)
	$(MPIFC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.f90 Makefile
	$(MPIFC) $(FFLAGS) -c $< -o $@

he3_kinds.o: ../src/he3_kinds.f90 Makefile
	$(MPIFC) $(FFLAGS) -c $< -o $@

legacy_anderson_mixing.o: ../src/legacy_anderson_mixing.f90 he3_kinds.o Makefile
	$(MPIFC) $(FFLAGS) -c $< -o $@

bulkgap.o: global_dec.o
init_calc.o: global_dec.o bulkgap.o
mpicalls.o: global_dec.o
interpol.o: global_dec.o
riccati.o: global_dec.o
getnewses.o: global_dec.o interpol.o riccati.o
packses.o: global_dec.o
iter_AA.o: global_dec.o getnewses.o packses.o legacy_anderson_mixing.o
iter.o: global_dec.o mpicalls.o getnewses.o packses.o iter_AA.o \
	iter_NN.f90 Iter_BR.f90 iter_BB.f90
qcv.o: global_dec.o init_calc.o iter.o

# Named mode targets clean first so objects built with different flag sets can
# never be mixed accidentally.
release:
	$(MAKE) clean
	$(MAKE) FFLAGS="$(COMMON_FFLAGS) $(RELEASE_FFLAGS) $(EXTRA_FFLAGS)" all

debug:
	$(MAKE) clean
	$(MAKE) FFLAGS="$(COMMON_FFLAGS) $(DEBUG_FFLAGS) $(EXTRA_FFLAGS)" all

show-config:
	@echo "MPI Fortran compiler: $(MPIFC)"
	@echo "Fortran flags:       $(FFLAGS)"
	@echo "Linker flags:        $(LDFLAGS)"
	@echo "Libraries:           $(LDLIBS)"

clean:
	rm -f $(OBJ) $(EXEC) *.mod

.PHONY: all release debug show-config clean
