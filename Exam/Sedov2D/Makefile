include ../sRules.make
CCFLAGS = $(OPT) $(specialrules)
LDFLAGS = -L. -lsedov2d  $(OPT)


VORO = ../../Voro
RKS = ../../RkS
COSMOS = ../../Cosmos
TIMER = ../../Timer
MPIAUX = ../../MpiAux
COLOR = ../../Colorize
OST = ../../OST
EXAM = ../
NR = ../../Cosmos/CosmosIC/


SEDOV2DLIB = libsedov2d.a

SEDOV2DOBJ = util.o sedov2d.o

INCLUDES = -I. -I../ -I../../  -I$(VORO) -I$(FFTW)/include -I$(OST) \
	-I$(RKS) -I$(COSMOS) -I$(TIMER) -I$(MPIAUX)  -I$(COLOR) \
	-I$(EXAM) -I$(NR)

LIBS = -L. -L$(VORO) -lmyvoro -L$(COLOR) -lcolor

default:
	$(MAKE) $(SEDOV2DLIB)


all:
	$(MAKE) $(SEDOV2DLIB)


new:
	$(MAKE) clean
	$(MAKE) $(SEDOV2DLIB)


$(SEDOV2DLIB): $(SEDOV2DOBJ)
	ar -rcv $@ $?
	ranlib $@


%.c:
	$(CC) $(CCFLAGS) -c $*.c


clean:
	-rm -f *.o *.a *.d core


.SUFFIXES: .exe .o .c .f .F .f90 .cu

.cc.o :
	$(CPP) $(CCFLAGS) $(INCLUDES) -c $<

.c.o :
	$(CC) $(CCFLAGS) $(INCLUDES) -c $<
