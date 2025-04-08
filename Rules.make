#
# Makefile for PM/Tree code pyul.exe - April 23, 2001 - JJD
#




RANLIB = ranlib




include ./sRules.make

#####################################
# System specific definitions
#####################################

#FFTWDIR = /user/kjhan/fftwfinal/
#FFTWDIR = /home/kjhan/local/
##
MKL_FFTWDIR = /home01/g053pcb/fftw_mkl_intel_openMPI_single/
MKL_LIB = /applic/compilers/intel/11.1/mkl/lib/em64t

#OPT = -DPGCC -mcmodel=medium -fast -fastsse
NVCC = nvcc --maxrregcount 128
#F90C = ifort
#FFTWDIR = /home/t075pcb/intel/
#OPT = -DINTEL -openmp -O2 -xSSE4.2 -m64

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
# -DWGROUPSIZE=nnn      - size of the subgroup for sequential I/O of data
# -DINCLUDE_TREE_FORCE - define if you want tree code corrections
# -DTREE               - define if you want tree code corrections
# -DDEBUG              - output debugging information
# -DBIT64              - use on machines using 64 bit addressing 
#                       - i.e. Compaq alpha
# -DAIX                - use on IBM SP3 machines
# -DOBSERVATION_MODE   - save halo data at every ObservationStride'th step
# -DSAVESLICE          - save slice density file at every time step
#############################################################


COOLLIBS = -L$(COOLINGDIR) -lget_lambda
SNFEEDBACKLIBS = -L$(SNFEEDBACKDIR) -lnew_ssp
EXTERNLIBS = -L$(EXTERNDIR) -lextforce
RTLIBS = -L$(RTDIR) -lraytay

#############################################################

INCLUDES = -I$(FFTW)/include -I$(FFT) -I$(MPIAUX) -I./ \
	-I$(RkS) -I$(Cosmos) -I$(CosmosIC) -I$(SPH) -I$(OSTDIR) \
	-I$(PARAMS) -I$(VORO) -I$(ParallelIO) -I$(KH) -I$(COLORLIB)  \
	-I$(RT) -I$(KEPLER) -I$(GLASS2D)

#####
#####
F90INCLUDES = -I$(CAMBDIR)/
COMFLAGS = $(OPT) $(specialrules)
FDFLAGS =  -DINCLUDE_TREE_FORCE  
CUFLAGS = $(CUOPT) -I$(CUDA)/include -L$(CUDA)/lib64/ -I$(CUDA_SDK)/common/inc -L$(CUDA_SDK)/lib 

# CUDA OPTIMIZATION
CUOPT = -O3  -use_fast_math  $(COMFLAGS)
# CUDA DIRECTORY
CUDA = /usr/local/cuda/
# CUDA SDK DIRECTORY
CUDA_SDK = /usr/local/cuda/C
CAMBLIBS = -L$(CAMBDIR) -lcamb #-pgf90libs
MPIAUXLIBS = -LMpiaux -lmpiaux


#LIBS = -L./FFT -lfft -L$(FFTW)/lib/ -lfftw3f_mpi  -lfftw3f -lfftw3f_threads   $(MPIAUXLIB) \
#	-L$(MPIAUX) -lmpiaux \
#	-L$(RMS) -lrms \
#	-L$(PARAMS) -lheader\
#	-L$(Cosmos) -lcosmology \
#	-L$(CosmosIC)  -ltwolpt \
#	-L$(CosmosRW) -lkjhrw \
#	-L$(RMS) -lrms \
#	-L$(TIMER) -ltimerutil\
#	-L$(SPH) -lsph \
#	-L$(VORO) -lmyvoro \
#	-L$(REC) -lrecfast \
#	-L$(OSTDIR) -lost  $(CAMBLIBS) \
#	-L$(DENRW) -ldenrw \
#	$(COOLLIBS) $(SNFEEDBACKLIBS) $(EXTERNLIBS) $(RTLIBS)\

LIBS = $(MPIAUXLIB) \
    -L$(PARAMS) -lheader\
    -L$(Cosmos) -lcosmology \
    -L$(CosmosIC)  -ltwolpt \
    -L$(CosmosRW) -lkjhrw \
    -L$(RkS) -lrks \
    -L$(TIMER) -ltimerutil\
    -L$(SPH) -lsph \
    -L$(VORO) -lmyvoro \
    -L$(REC) -lrecfast \
    -L$(DENRW) -ldenrw \
    -L$(MPIAUX) -lmpiaux \
    -L$(ParallelIO) -lparallelIO \
    -L$(KH) -lkh \
    -L$(VORO) -lmyvoro \
	-L$(COLORLIB) -lcolor\
	-L$(RT) -lrt \
	-L$(KEPLER) -lkp \
	-L$(GLASS2D) -lglass2d \
	-L$(EXAM) -lexam \
    -L$(OSTDIR) -lost \
    $(COOLLIBS) $(SNFEEDBACKLIBS) $(EXTERNLIBS) $(RTLIBS)\
    -L./FFT -lfft -L$(FFTW)/lib/ -lfftw3f_mpi  -lfftw3f -lfftw3f_threads   \
	$(CAMBLIBS) 


		# -L$(CUDA)/lib64 -lcudart
#	$(MKL_LIB)/libmkl_intel_lp64.a -Wl,--start-group -L$(MKL_LIB) -lmkl_cdft_core -lmkl_blacs_openmpi_lp64\
#	-lmkl_sequential  -lmkl_core -lmkl_intel_lp64 -Wl,--end-group





#CDFLAGS = -DWGROUPSIZE=500 -DNMEG=21500L -DINCLUDE_TREE_FORCE \

CDFLAGS = -DNMEG=8000L -DINCLUDE_TREE_FORCE \
        -D_LARGE_FILES -DSAVESLICE  -DPMSEEDFORCE   -DUSE_GPU\
		#-DTSC_OLD -DOLD_FDA4


FFLAGS = $(FDFLAGS) $(COMFLAGS)  
CFLAGS = $(CDFLAGS)  $(COMFLAGS) 
#LDFLAGS = $(OPT)  -nofor-main $(specialrules)
LDFLAGS = $(OPT)  $(specialrules) 
CAMBLIBS = -L$(CAMBDIR) -lcamb



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
.SUFFIXES: .exe .o .c .f .F .f90 .cu

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
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS)

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
