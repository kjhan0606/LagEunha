#ifndef GFS_MATERIAL_TRACE_H
#define GFS_MATERIAL_TRACE_H
#include <math.h>
#ifdef __CUDACC__
#define GFS_TRACE_HD __host__ __device__
#else
#define GFS_TRACE_HD
#endif
#define GFS_MATERIAL_TRACE_N 4
/* Fixed material-reference face space. A coefficient multiplier mu represents
 * actual reference traction t=sum phi_k eta_k, G eta=mu, G=int phi phi^T dS0.
 * Its load is dual to the SAME kinematic trace. No torque/heat correction.
 * At a deformed face t_physical=t_reference/(dS/dS0). */
GFS_TRACE_HD static inline int gfs_material_trace_factor(int n,const double G[4][4],double L[4][4]){
    if(n<1 || n>4)return 0;
    for(int i=0;i<n;i++)for(int j=0;j<=i;j++){
        double a=G[i][j];for(int k=0;k<j;k++)a-=L[i][k]*L[j][k];
        if(i==j){if(!(a>0) || !isfinite(a))return 0;L[i][j]=sqrt(a);}else L[i][j]=a/L[j][j];
    }return 1;
}
GFS_TRACE_HD static inline void gfs_material_trace_solve(int n,const double L[4][4],double rhs[4]){
    for(int i=0;i<n;i++){for(int j=0;j<i;j++)rhs[i]-=L[i][j]*rhs[j];rhs[i]/=L[i][i];}
    for(int i=n-1;i>=0;i--){for(int j=i+1;j<n;j++)rhs[i]-=L[j][i]*rhs[j];rhs[i]/=L[i][i];}
}
#undef GFS_TRACE_HD
#endif
