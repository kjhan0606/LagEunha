#########################################################
####      KIAS INTEL Compiler                       #####
#########################################################
AR = ar rcv
### intel mpi compilers
#FC = mpiifort -nofor-main 
#CC = mpiicc -cc=icx
#CPP = icx
### openmpi compilers
FC = mpifort -nofor-main 
CC = mpicc  -w2
CPP = icc -w2


F90C = ifort -nofor-main 
### sgi compilers
#FC = ftn
#CC = cc
#F90C = ftn


##########
# optimization options
#OPT = -O3     -qopenmp #-xavx #-E
OPT = -g       -qopenmp #-xavx #-E

##########################################################################
### (1) This is for PM-type N-body
#S_TYPE = -DGOTPM -DPMonly

### (2) This is for GOTPM-type N-body
S_TYPE = -DGOTPM

### (3) This is for Hydro SPH simulation which is, however, not yet implemented.
#S_TYPE =  -DSPH


### (3) This is for Hydro VSPH simulation which is, however, not yet implemented.
#S_TYPE =  -DVSPH
##########################################################################

### THIS IS THE TYPICAL GOTPM COMPILING OPTION ###
#specialrules = -DINTEL   -DGOTPM   -DXYZDBL -DDEBUG
#specialrules = -DINTEL   -DGOTPM   -DDEBUG

specialrules = -DINTEL   -DXYZDBL  $(S_TYPE)    -DDEBUG   #  -Wmissing-prototypes #-DGOTPM # -DDEBUG

FFTW = /home/kjhan/local

#FFTW = ./fftw/
