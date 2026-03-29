//*******************
// The number of interation to determine 
// the boundedness of member particles
#define BOUNDITER 4
//*******************
//


//*******************
// fractional error to measure 
// the core density in ratio
#define COREDENRESOLUTION (1.e-3) 
//*******************


//*******************
// lowest stellar-density peak
// in unit of h^2 Msun/ckpc^3
#define PEAKTHRESHOLD 2000.L
//*******************

//--------------------------------
// please tune these six parameters 
// for higher resolution simulations 
//
//*******************
// the minimum number of star/dm 
// particles to identify a core 
#define MINCORENMEM 30 
//*******************

//*******************
// the number of neighbors to build 
// the neighbor network 
// the bigger the better
// it should be no larger than MAX_NUM_NEAR in tree.h
#define NUMNEIGHBOR 32 
//*******************

//*******************
// the number of iso-den division of non-core 
// particles 
#define NSHELLDIVIDE 10 
//*******************

//*******************
// the separation limit of peaks 
// to merge in cMpc/h 
#define MERGINGPEAKLENGTH 4.e-3  
//*******************

//*******************
// the minimum stellar mass for a core 
#define MINCORESTARMASS -1 
//*******************

//*******************
// the minimun stellar mass of the FoF halo 
// for galaxy finding with stellar density 
#define MINSTELLARMASS 2.e6  
//*******************
//--------------------------------

//*******************
// Maximum number of cores 
#define MAXNUMCORE 1000000
//*******************


//*******************
// (obsolete) The number of nearby stars 
// to measure stellar density
#define NUMNEARDEN 10
//*******************

//*******************
// (obsolete) The number of stellar neighbors
// for the core detection to find density peaks 
#define NUMSTELLARNEIGHBORS 30 
//*******************


//*******************
// (obsolete) minimum smoothing length 
// for stellar density in cMpc/h
#define MIN_CONST_R_SMOOTHING 0.010
//*******************

//*******************
// The cellsize for TSC of the stellar density 
// in unit of cMpc/h
#define TSC_CELL_SIZE 0.004
//*******************

//*******************
// Gaussian Smoothing Length 
// for stellar density in unit of cMpc/h. 
// This should be not smaller 
// than 2*TSC_CELL_SIZE.
// Note that the smoothing may 
// lower the stellar density 
// and PEAKTHRESHOLD should be 
// lowered accordingly.
#define Gaussian_Smoothing_Length 0.008
//*******************
//
//*******************
// mininum density
#define DENFLOOR 1.
//*******************
//
//
//
// the (even) number of cell boundary buffer. 
#define NCELLBUFF 10
//*******************
//
//
//
//
//
// the size of deep linked lists (omp parallelized)
// if the number of particles linked to a cell is smaller
// than DEEPSIZE, the omp parallelization is located outside the loop.
// In the other cases, the omp parallelization goes down the loop.
#define DEEPSIZE 1024
//*******************
//
//
// The lower Number of particles for the OMP-parallized FoF
#define NOMPFoF 500000
//*******************
//
//
//
// the linking length to finalize the membership
#define FOFLINK4MEMBERSHIP 0.005
//*******************
//
// the maximum number of linking in water shedding to find core density
// It should be sufficiently larger than 
// the maximum number of core particles.
#define MAXNUMWATERSHEDDING 100000000L
//*******************
//
//
// the maximum number of threads
#define MAXTHREADS 64
