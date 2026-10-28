#ifndef LAG_MOMENT_GRAVITY_QUADRATURE_H
#define LAG_MOMENT_GRAVITY_QUADRATURE_H
#include <math.h>
#include <stddef.h>
#ifdef __CUDACC__
#define MGS_HD __host__ __device__
#else
#define MGS_HD
#endif
/* Samples: x,y,dm,vx,vy,ax,ay. Output: force2,D4,Pi(xx,xy,yy),work.
 * Acceleration-only source; no position or potential-energy kick.
 * Invalid inputs/overflow leave output unchanged. */
static MGS_HD int mgs_integrate(const double *samples,size_t count,double *out){
    if(!samples || !out || count==0)return 1;
    double s[10]={0};
    for(size_t i=0;i<count;i++){
        const double *r=samples+7*i;
        for(int k=0;k<7;k++)if(!isfinite(r[k]))return 1;
        if(!(r[2]>0))return 1;
        double fx=r[2]*r[5],fy=r[2]*r[6];
        s[0]+=fx;s[1]+=fy;
        s[2]+=r[0]*fx;s[3]+=r[0]*fy;s[4]+=r[1]*fx;s[5]+=r[1]*fy;
        s[6]+=2*r[3]*fx;s[7]+=r[3]*fy+r[4]*fx;s[8]+=2*r[4]*fy;
        s[9]+=r[3]*fx+r[4]*fy;
    }
    for(int k=0;k<10;k++)if(!isfinite(s[k]))return 2;
    for(int k=0;k<10;k++)out[k]=s[k];
    return 0;
}
#undef MGS_HD
#endif
