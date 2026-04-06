#
# Makefile for PM/Tree code pyul.exe - April 23, 2001 - JJD
#

#--- Executable

PROGS =  eunha2.exe 

#--- Objects to build

OBJS		= eunha2.o 

OS = LINUX_GCC



#--- Targets

CAMBDIR = ../Cosmos/CosmosIC/CAMB_Sources
MPIAUX = ../MpiAux
PARAMS = ../Params
Cosmos = ../Cosmos
CosmosIC = ../Cosmos/CosmosIC
CosmosRW = ../Cosmos/CosmosRW
SPH = ../SPH
RkS = ../RkS
OSTDIR = ../OST
REC = ../RECFAST
COOLINGDIR = ../SPH/JShin/cooling
SNFEEDBACKDIR = ../SPH/JShin/ssp
EXTERNDIR = ../SPH/JShin/Extern
TIMER = ../Timer
DFFT = ../FFT
DENRW = ../DenRW
ParallelIO = ../ParallelIO
EXAM = ../Exam
KH = ../Exam/KH
VORO = ../Voro
COLORLIB = ../Colorize
RT = ../Exam/RT
KEPLER = ../Exam/Kepler
GLASS2D = ../Exam/MkGlass2D


SUBDIRS = $(MPIAUX) $(CAMBDIR) $(PARAMS) $(Cosmos) $(CosmosIC) \
		$(SPH) $(VORO) $(RkS) $(OSTDIR) $(REC) $(COOLINGDIR) \
		$(SNFEEDBACKDIR) $(EXTERNDIR) $(RTDIR) $(CosmosRW) \
		$(TIMER) $(DFFT) $(DENRW) $(ParallelIO) $(KH) $(VORO)\
		$(COLORLIB)  $(RT) $(KEPLER) $(GLASS2D) $(EXAM)

all:
	$(MAKE) misc

this: $(OBJS)
	rm -f $(PROGS) 
	$(MAKE) $(PROGS)
new: 
	$(MAKE) clean
	$(MAKE) all

clean: 
	rm -f $(OBJS)
	rm -f $(OBJPOW)
	rm -f $(OBJVIEW)
	rm -f $(PROGS)
	rm -f $(MKPOW)
	rm -f $(VIEWPOW)
	rm -rf $(MKPOW) $(VIEWPOW) $(MKPFILE) $(EXAMREAD)

MKPOW = mkpower.exe
OBJPOW = mkpower.o
VIEWPOW = viewpower.exe
OBJVIEW = viewpower.o
MKPFILE = mkpfile.exe
MKTRACK = mktrack.exe
GPUTEST = gputest.exe
EXAMREAD = sync2txt.exe


mkview:
	$(CC) -c -DNMEG=1000L mkviewer.mod2.c mktrack.c -g
	$(F90C) -o $(MKTRACK) mktrack.o mkviewer.mod2.o ratint.f Memory2.o -g
	

misc:
	$(F90C) -c viewpower.f90 -I./$(CAMBDIR) $(OPT) $(specialrules)
	$(F90C) -o $(VIEWPOW) viewpower.o -I./$(CAMBDIR) -L./$(CAMBDIR) -lcamb $(OPT) $(specialrules)
	$(F90C) -c mkpower.f90 -I./$(CAMBDIR) $(OPT) $(specialrules)
	$(CC) -c getsimparm.c $(OPT) $(specialrules) -lm -I$(FFTW)/include -I$(RkS)  -I$(PARAMS)
	$(FC) -o $(MKPOW) mkpower.o getsimparm.o -L$(PARAMS) -lheader -L$(MPIAUX) -lmpiaux -I./$(CAMBDIR) -L./$(CAMBDIR) -lcamb $(OPT) $(specialrules) -L$(Cosmos) -lcosmology -L$(CosmosIC) -ltwolpt
	$(CC) -o $(MKPFILE) mkpfile.c -I$(FFTWDIR)/include  -L$(PARAMS) -lheader  -L$(MPIAUX) -lmpiaux  $(OPT) $(specialrules) -lm   -I$(RkS) -I$(PARAMS) -I$(FFTW)/include -L$(Cosmos) -lcosmology -L$(CosmosIC) -ltwolpt -L$(DENRW) -ldenrw # -l$(VORO) -lmyvoro

# sync2txt and econ require eunha2.o (my_malloc/my_free) — build from top-level src/
heavy:
	$(CC) -c sync2txt.c  $(OPT) $(specialrules) -lm $(INCLUDES) $(MPI_LINKS)
	$(FC) -o $(EXAMREAD) sync2txt.o ../eunha2.o -qopenmp $(MPI_LINKS) $(INCLUDES)  $(LDFLAGS) $(LIBS) -lm
	$(CC) -c econ.c  $(OPT) $(specialrules) -lm $(INCLUDES) $(MPI_LINKS)
	$(FC) -o econ.exe econ.o ../eunha2.o -qopenmp $(MPI_LINKS) $(INCLUDES)  $(LDFLAGS) $(LIBS) -lm


gputest:
	rm -rf {GPUTEST}
	gcc -c gputest.c   $(OPT) $(specialrules)
	gcc -c Treewalk.c   $(OPT) $(specialrules)
	$(NVCC) -o $(GPUTEST) force_spline.o gputest.o Memory2.o gpu.Treewalk.mod2.cu  Treewalk.o \
		-lm -L$(CUDA)/lib64/ -L$(CUDA_SDK)/lib -I$(CUDA)/include \
		-I$(CUDA_SDK)/common/inc  timerutil.o  $(OPT) $(specialrules)   #-deviceemu

default: all

#$(MKPOW): $(OBJPOW)



include Rules.make
include sRules.make
