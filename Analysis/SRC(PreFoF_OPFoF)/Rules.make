#
# Makefile for PM/Tree code pyul.exe - April 23, 2001 - JJD
#


RANLIB = ranlib


include ../../sRules.make



#####################################
# System specific definitions
#####################################

#########################################################
####      KIAS PG Compiler                          #####
#########################################################
#AR = ar rcv
#FC = mpiifort
#CC = mpiicc
#F90C = icc
#FFTWDIR = /home01/e994pcb/local/
#OPT = -DPGCC -mcmodel=medium -tp nehalem-64 -fast -mp -fastsse
#OPT = -g -DINTEL -mcmodel=medium 

#########################################################
####      KISTI IBM SP MACHINES                     #####
#########################################################
#AR = ar -X64 rcv
#FC = mpxlf_r
#CC = mpxlc_r
#F90C = mpxlf_r
#FFTWDIR = /user/kjhan/fftwfinal/
#OPT = -O3  -q64 -qtune=pwr5 -qarch=pwr5

#########################################################
####      KISTI IBM SP MACHINES                     #####
#########################################################
#AR = ar rcv
#FC = mpif77
#CC = mpicc
#F90C = ifort
#FFTWDIR = /user/kjhan/fftwfinal/
#OPT = -O1 -xP -ip


#############################################
# List of compilation directives
#############################################
# 
# -DNMEG=nnn      - number of megabytes of storage per processor - default=40
# -DWGSIZE=nnn      - size of the subgroup for sequential I/O of data
# -DINCLUDE_TREE_FORCE - define if you want tree code corrections
# -DTREE               - define if you want tree code corrections
# -DDEBUG              - output debugging information
# -DBIT64              - use on machines using 64 bit addressing 
#                       - i.e. Compaq alpha
# -DAIX                - use on IBM SP3 machines
# -DOBSERVATION_MODE   - save halo data at every ObservationStride'th step
# -DSAVESLICE          - save slice density file at every time step
#############################################################

FFTWDIR = /home/kjhan/local/

PARAMS = ../../Params/
MPIAUX = ../../MpiAux/
COSMOS = ../../Cosmos/
TWOLPT = ../../Cosmos/CosmosIC/
RMS = ../../RMS
BASE = ../../


INCLUDES = -I$(FFTWDIR)/include

#LIBS = -L$(FFTWDIR)/lib/ -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw   $(CAMBLIBS) 
#####
#####
F90INCLUDES = -I./$(CAMBDIR)
COMFLAGS = -DINDEX  -DOLD  $(OPT) $(specialrules)


LIBS = -L$(PARAMS) -lheader -L$(MPIAUX) -lmpiaux -L$(COSMOS) -lcosmology -L$(TWOLPT) -ltwolpt 

FDFLAGS =  -DINCLUDE_TREE_FORCE  

CDFLAGS = -DWGSIZE=20 -DNMEG=3000L -DINCLUDE_TREE_FORCE -I./$(PARAMS)  -I$(PARAMS) -I$(MPIAUX) -I$(COSMOS) -I$(TWOLPT) -I$(BASE)\
        -I$(RMS) -D_LARGE_FILES -DSAVESLICE  -DPMSEEDFORCE # -DUSE_MASTER		#-DTSC_OLD -DOLD_FDA4


FFLAGS = $(FDFLAGS) $(OPT) $(COMFLAGS)  
CFLAGS = $(OPT) $(CDFLAGS)  $(COMFLAGS) 
LDFLAGS = $(OPT)  $(specialrules)
CAMBLIBS =



#################################
# Compaq Alpha
#################################
#FC = f77
#CC = cc 
#INCLUDES = -I/home/kjhan/dolphin/fftw/include
#FDFLAGS = -DBIT64 -DINCLUDE_TREE_FORCE -DTREE -DCOMPACT
#DFLAGS = -DBIT64 -DNMEG=256 -DINCLUDE_TREE_FORCE -DTREE -DCOMPACT -DTREEFIX
#FFLAGS = $(FDFLAGS)  -fast -nofor_main
#CFLAGS = $(DFLAGS)  -fast
#LDFLAGS =  -fast -nofor_main
#LIBS = -L/home/kjhan/dolphin/fftw/lib -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw -lm -lz -lmpi -lm
##################################

#--- C Compiler information
#  Leave the rest untouched


#--- Suffix-based compilation rules
.SUFFIXES: .exe .o .c .f .F .f90

#rules to build binary from source


.c.o :
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

.f90.o :
	$(F90C) $(FFLAGS) $(INCLUDES) $(F90INCLUDES) -c $<

.f.o :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.for.o :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.F.o :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.cu.o :
	$(NVCC) $(CUFLAGS) $(INCLUDES) -c $<

.c.exe :
	$(CC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS)

.f.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS) $(INCLUDES)

.f90.exe :
	$(F90C) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS) $(INCLUDES)

.for.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS) $(INCLUDES)

.F.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS) $(INCLUDES)

.cu.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS)

#--- Targets
