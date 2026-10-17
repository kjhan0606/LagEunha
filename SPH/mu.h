#define nDelta 128
#define nT 128
#define nTMU 128

#ifdef _MAIN_MU 
#define EXTERN
#else
#define EXTERN extern
#endif


EXTERN double aMU[nT*nDelta];
EXTERN double aTMU[nT*nDelta];
EXTERN double aT[nTMU*nDelta];


#define aMU(i,j) aMU[(i)*nDelta+(j)]
#define aTMU(i,j) aTMU[(i)*nDelta+(j)]
#define aT(i,j) aT[(i)*nDelta+(j)]
