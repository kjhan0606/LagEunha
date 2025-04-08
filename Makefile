#
# Makefile for PM/Tree code pyul.exe - April 23, 2001 - JJD
#

#--- Executable

PROGS =  eunha2.exe 

#--- Objects to build

OBJS		= eunha2.o 
#			observermain.o \
#			observerlightcone.obs0.o observerlightcone.obs1.o observerlightcone.obs2.o \
#			observerlightcone.obs3.o observerlightcone.obs4.o observerlightcone.obs5.o \
#			observerlightcone.obs6.o observerlightcone.obs7.o\
#			PencilBeamlightcone.o observerlightcone.obsCen.o\
			# fda4_new.o fda4_split.o tsc_c.o

#			fmod.o fda4_new.o fda4_split.o tsc_c.o \
#indexing.o PMPI_Sendrecv_Replace.o  #kjhtpm.mod2.o
#hypersort.o\
#slicelightcone.final.o
#			kjhtpm.final.o indexing.o PMPI_Sendrecv_Replace.o  #kjhtpm.mod2.o
#			findindex.o maxextent.o glacial.o survey.o surveywrite2.o
#			memory.o	tpm.o domdecomp.o  kjhgetvisitors.o \
# 			tree.o  quickbin.o group.o groupgrav.o  \
#			kjhtpm.o Treewalk.mod6.o force_spline.mod2.o
#	hypersort.o \

OS = LINUX_GCC



#--- Targets

#CAMBDIR = ./CAMB
CAMBDIR = ./Cosmos/CosmosIC/CAMB_Sources
MPIAUX = ./MpiAux
PARAMS = ./Params
Cosmos = ./Cosmos
CosmosIC = ./Cosmos/CosmosIC
CosmosRW = ./Cosmos/CosmosRW
SPH = ./SPH
RkS = ./RkS
OSTDIR = ./OST
REC = ./RECFAST
COOLINGDIR = ./SPH/JShin/cooling
SNFEEDBACKDIR = ./SPH/JShin/ssp
EXTERNDIR = ./SPH/JShin/Extern
TIMER = ./Timer
FFT = ./FFT
DENRW = ./DenRW
ParallelIO = ./ParallelIO
EXAM = ./Exam
KH = ./Exam/KH
VORO = ./Voro
COLORLIB = ./Colorize
RT = ./Exam/RT
KEPLER = ./Exam/Kepler
GLASS2D = ./Exam/MkGlass2D


SUBDIRS = $(MPIAUX) $(CAMBDIR) $(PARAMS) $(Cosmos) $(CosmosIC) \
		$(SPH) $(VORO) $(RkS) $(OSTDIR) $(REC) $(COOLINGDIR) \
		$(SNFEEDBACKDIR) $(EXTERNDIR) $(RTDIR) $(CosmosRW) \
		$(TIMER) $(FFT) $(DENRW) $(ParallelIO) $(KH) $(VORO)\
		$(COLORLIB)  $(RT) $(KEPLER) $(GLASS2D) $(EXAM)

all: 
	for subdir in $(SUBDIRS); do \
		(cd $$subdir && $(MAKE) );\
	done 
	$(MAKE) this

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
	for subdir in $(SUBDIRS); do \
		(cd $$subdir && $(MAKE) clean);\
	done 

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
	rm -f $(VIEWPOW) 
	rm -f $(OBJVIEW) 
	rm -f $(MKPFILE) 
	rm -f $(MKTRACK) 
	$(F90C) -c viewpower.f90 -I./$(CAMBDIR) $(OPT) $(specialrules)
	$(F90C) -o $(VIEWPOW) viewpower.o -I./$(CAMBDIR) -L./$(CAMBDIR) -lcamb $(OPT) $(specialrules)
	$(F90C) -c mkpower.f90 -I./$(CAMBDIR) $(OPT) $(specialrules)
	$(CC) -c getsimparm.c $(OPT) $(specialrules) -lm -I$(FFTW)/include -I$(RkS)  -I$(PARAMS)
	$(FC) -o $(MKPOW) mkpower.o getsimparm.o Params/header.o -L$(MPIAUX) -lmpiaux -I./$(CAMBDIR) -L./$(CAMBDIR) -lcamb $(OPT) $(specialrules) -L$(Cosmos) -lcosmology -L$(CosmosIC) -ltwolpt
	$(CC) -o $(MKPFILE) mkpfile.c -I$(FFTWDIR)/include  -L$(PARAMS) -lheader  -L$(MPIAUX) -lmpiaux  $(OPT) $(specialrules) -lm   -I$(RkS) -I$(PARAMS) -I$(FFTW)/include -L$(Cosmos) -lcosmology -L$(CosmosIC) -ltwolpt -L$(DENRW) -ldenrw # -l$(VORO) -lmyvoro
	$(CC) -c sync2txt.c  $(OPT) $(specialrules) -lm $(INCLUDES) $(MPI_LINKS)  $(LDFLAGS) $(LIBS) 
	$(FC) -o $(EXAMREAD) sync2txt.o -qopenmp $(MPI_LINKS) $(INCLUDES)  $(LDFLAGS) $(LIBS) -lm

#	$(MAKE) mkview
#	$(CC) -o readtest readtest.c Memory2.o timerutil.o  header.o   $(OPT) $(specialrules)


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
