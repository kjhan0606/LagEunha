include ../sRules.make
CCFLAGS = $(OPT) $(specialrules)

VORO = ../../Voro
RKS = ../../RkS
COSMOS = ../../Cosmos
TIMER = ../../Timer
MPIAUX = ../../MpiAux
COLOR = ../../Colorize
OST = ../../OST
EXAM = ../

CYLLIB       = libcyl.a

CYLOBJ       = util.o cylinder.o

INCLUDES = -I. -I../ -I../../  -I$(VORO) -I$(FFTW)/include -I$(OST)\
		-I$(RKS) -I$(COSMOS) -I$(TIMER) -I$(MPIAUX) -I$(COLOR)  \
		-I$(EXAM)

default:
	$(MAKE) $(CYLLIB)

util.o: util.c cylinder.h ../gfs_material_radial_image.h
cylinder.o: cylinder.c cylinder.h material_boundary.inc material_boundary.h ../gfs_material_open_stream.h

$(CYLLIB): $(CYLOBJ)
	ar -rcv $@ $?
	ranlib $@

clean:
	-rm -f *.o *.a *.d core

.SUFFIXES: .exe .o .c

.c.o :
	$(CC) $(CCFLAGS) $(INCLUDES) -c $<
