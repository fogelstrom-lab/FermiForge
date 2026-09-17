FFLAGS  = -fdefault-real-8 -O4 -fallow-argument-mismatch -march=native -ftree-vectorizer-verbose=2

FCODES  = iter.f gap.f riccati.f linpacks.f trajectories.f init_calc.f mpicalls.f interpol.f
OBJS    = iter.o gap.o riccati.o linpacks.o trajectories.o init_calc.o mpicalls.o interpol.o
FC      = gfortran
MPIcomp	= mpifort

qcv      :  qcv.o  $(OBJS) $(FCODES) Makefile
	$(MPIcomp) $(FFLAGS) -o qcv qcv.o $(OBJS)
clean	:
	rm *.o qcv 

.f.o    : Makefile $(FCODES)
	$(FC) $(FFLAGS) -c $*.f

iter.o        : iter.f iter.dat
qcv.o         : qcv.f qcv.dat
init_calc.o   : init_calc.f qcv.dat
gap.o         : gap.f gap.dat
riccati.o     : riccati.f qcv.dat
gap.o         : gap.f gap.dat
linpacks.o    : linpacks.f

