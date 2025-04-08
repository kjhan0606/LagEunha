#CAMB sources makefile
#Edit for your compiler

include ../../../sRules.make

FFLAGS = $(OPT) $(specialrules)

#Intel ifort , -openmp toggles mutli-processor:
#Beware of openmp bugs in version 10.0 (10.1 seems OK)
#F90C     = ifort
#FFLAGS = -ip -O2 -vec_report0 -W0 -WB -fpp2
#include ../../../Rules.make
#FFLAGS = -openmp -O2 -ip -W0 -WB -fpp2 -vec_report0
#FFLAGS = -openmp -O2 -ip -W0 -WB -fpp2 -vec_report0
#F90C = pgf90
#FCFLAGS = -O3 -mcmodel=medium
#CCFLAGS = -O3 -mcmodel=medium -I../../../ 

#CC = mpicc
#FC = mpif77

#Sun, single processor:
#F90C     = f90
#FFLAGS = -O2

# Intel 9 on IA-64 (eg. COSMOS)
# (do "module load icomp90" before compiling)
#F90C = ifort
#FFLAGS = -openmp -fpp2 -w -O3 -ip -mP2OPT_hlo_prefetch=F

#Intel ifc, add -openmp for multi-processor (some have bugs):
#F90C     = ifc
#FFLAGS = -O2 -Vaxlib -ip -W0 -WB -quiet -fpp2
#some systems can can also add e.g. -tpp7 -xW

#SGI, -mp toggles multi-processor. Use -O2 if -Ofast gives problems.
#F90C     = f90
#FFLAGS  = -Ofast -mp

#Digital/Compaq fortran, -omp toggles multi-processor
#F90C    = f90
#FFLAGS  = -omp -O4 -arch host -math_library fast -tune host -fpe1

#Absoft ProFortran, single processor:
#F90C     = f95
#FFLAGS = -O2 -cpu:athlon -s -lU77 -w -YEXT_NAMES="LCS" -YEXT_SFX="_"

#NAGF95, single processor:
#F90C     = f95
#FFLAGS = -DNAGF95 -O3 -mismatch

#Gfortran compiler
#F90C = gfc
#FFLAGS = -O2 

#G95 compiler
#F90C   = g95
#FFLAGS = -O2

#Sun V880
#F90C = mpf90
#FFLAGS =  -O4 -openmp -ftrap=%none -dalign -DMPI

#Sun parallel enterprise:
#F90C     = f95
#FFLAGS =  -O2 -xarch=native64 -openmp -ftrap=%none
#try removing -openmp if get bus errors. -03, -04 etc are dodgy.

#IBM XL Fortran, multi-processor (run gmake)
#F90C     = xlf90_r
#FFLAGS  = -qsmp=omp -qsuffix=f=f90:cpp=F90 -O3 -qstrict -qarch=pwr3 -qtune=pwr3




FEEDBACKLIB       = libextforce.a

#Shouldn't need to change anything else...

F90FLAGS      = $(FFLAGS)
#FC            = $(F90C)

FEEDBACKOBJ       = externalforce.o readextforce.o jhsrw.o


CFLAGS = -I../../..  -I../../../RkS  -I$(FFTW)/include $(OPT) $(specialrules)

default: $(FEEDBACKLIB)

all: $(FEEDBACKLIB)


$(FEEDBACKLIB): $(FEEDBACKOBJ)
	ar -r $@ $?

#CC = mpicc

#FC = ifort

#	ranlib $(COOLLIB)

%.o: %.for
	$(FC) $(INCLUDES) $(FFLAGS) -c $*.for

%.o: %.c
	$(CC) $(INCLUDES) $(CFLAGS) -c $*.c

clean:
	-rm -f *.o *.a *.d core *.mod

new:
	$(MAKE) clean
	$(MAKE) all
	 
