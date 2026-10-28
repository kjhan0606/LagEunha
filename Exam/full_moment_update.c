#include "full_moment_update.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int full_moment_build_candidate(int n,double h,const double *base,const double *mass,
 const double *stage,const double *gen,const double *z,const long double *Psi,const double *Psi_rate,
 const double *factor,double *out,long double *out_Psi,double *out_L){
 if(n<1||!isfinite(h)||h<=0.||!base||!mass||!stage||!gen||!z||!Psi||!Psi_rate||!factor||!out||!out_Psi||!out_L)return 1;
 double *state=malloc((size_t)n*20*sizeof(double)),*L=malloc((size_t)n*4*sizeof(double));
 long double *potential=malloc((size_t)n*sizeof(long double));int status=0;
 if(!state||!L||!potential){free(state);free(L);free(potential);return 2;}
 for(int i=0;i<n&&!status;i++){
  const double *a=base+20*i,*f=factor+63*i,*r=stage+94*i;
  if(!isfinite(mass[i])||mass[i]<=0.||!isfinite(Psi[i])||!isfinite(Psi_rate[i])||r[93]!=0.||f[62]!=0.){status=3;break;}
  for(int j=0;j<63;j++)if(!isfinite(f[j]))status=3;
  if(!(f[0]>0.&&f[3]>0.&&f[1]==0.&&f[51]>0.&&f[53]<128*DBL_EPSILON)){status=5;break;}
  if(memcmp(f+4,mass+i,sizeof(double))||memcmp(f+18,f+19,sizeof(double))){status=4;break;}
  for(int k=1;k<18;k++)if(k<13||k==17){
   double expected=a[k+2]+h*r[k];
   if(memcmp(f+4+k,&expected,sizeof(double))){status=4;break;}
  }
  double *b=state+20*i;b[0]=a[0]+h*gen[2*i];b[1]=a[1]+h*gen[2*i+1];b[2]=a[2]+h*z[i];
  for(int k=0;k<17;k++)b[k+3]=f[k+5];
  potential[i]=Psi[i]+(long double)h*(long double)Psi_rate[i];
  for(int k=0;k<4;k++)L[4*i+k]=f[k];
  for(int k=0;k<20;k++)if(!isfinite(b[k]))status=3;
  if(!isfinite(potential[i]))status=3;
 }
 if(!status){memcpy(out,state,(size_t)n*20*sizeof(double));memcpy(out_Psi,potential,(size_t)n*sizeof(long double));memcpy(out_L,L,(size_t)n*4*sizeof(double));}
 free(state);free(L);free(potential);return status;
}
