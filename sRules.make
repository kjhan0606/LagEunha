#########################################################
####      KIAS INTEL Compiler                       #####
#########################################################
AR = ar rcv
### intel mpi compilers
#FC = mpiifort -nofor-main 
#CC = mpiicc -cc=icx
#CPP = icx
### openmpi compilers
# module intel/compiler/2021.3.0 openmpi/4.1.6_INTEL intel/compiler-rt/2021.3.0
#FC = mpifort -nofor-main 
#CC = mpicc  -w2 -debug
#CPP = icc -w2 -debug
#F90C = ifort -nofor-main 
#OPT = -g -O0 -qopenmp #-xavx #-E

# 2025/11/13 grammar version
FC = mpiifx -nofor-main 
CC = mpiicx   -debug 
CPP = icx -debug
F90C = ifx -nofor-main 
OPT = -g -O0 -qopenmp #-xavx #-E

### openmpi compilers


### sgi compilers
#FC = ftn
#CC = cc
#F90C = ftn


##########
# optimization options
#OPT = -O3     -qopenmp #-xavx #-E

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
