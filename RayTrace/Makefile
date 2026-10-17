#CAMB Makefile
#Edit for your compiler
#Note there are many ifc versions, some of which behave oddly

include ../sRules.make
CCFLAGS = $(OPT) $(specialrules)
LDFLAGS = -L. -lmyvoro


VOROLIB = ../Voro
VOROINC = ../Voro

TIMERLIB       = libmyrt.a

TIMEROBJ       = raytrace.o  hpxvertex.o


INCLUDES = -I. -I../ -I$(VOROINC) -I$(FFTW)/include -I../RkS/


LIBS = -L. -L$(VOROLIB)

default: 
	$(MAKE) $(TIMERLIB)


all: 
	$(MAKE) $(TIMERLIB)

new: 
	$(MAKE) clean  
	$(MAKE) $(TIMERLIB)



$(TIMERLIB): $(TIMEROBJ)
	ar -rcv $@ $?


%.c: 
	$(CC) $(CCFLAGS) -c $*.c

%.f90: 
	$(FC) $(FCFLAGS) -c $*.f90

clean:
	-rm -f *.o *.a *.d core 


#--- Suffix-based compilation rules
.SUFFIXES: .exe .o .c .f .F .f90 .cu

#rules to build binary from source

.cc.o : 
	$(CPP) $(CCFLAGS) $(INCLUDES) -c $<

.c.o : 
	$(CC) $(CCFLAGS) $(INCLUDES) -c $<

.c.exe :
	$(CC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(TESTOBJS) $(LIBS)

.f90.o : 
	$(F90C) $(FCFLAGS) $(INCLUDES) -c $<

.f90.exe :
	$(F90C) $(LDFLAGS) $(INCLUDES) -o $*.exe $(TESTOBJS) $(LIBS)


#include ../sRules.make
