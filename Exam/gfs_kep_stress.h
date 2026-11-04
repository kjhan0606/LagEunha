#ifndef GFS_KEP_STRESS_H
#define GFS_KEP_STRESS_H
#include <math.h>
#ifdef __CUDACC__
#define GFS_STRESS_HD static __host__ __device__ inline
#else
#define GFS_STRESS_HD static inline
#endif
/* The energies use the ORIGINAL unlimited centroid LSQ objective weights,
 * A/|r|^(d+1). No fitted threshold, extra moment or physical viscosity. */
GFS_STRESS_HD double gfs_kep_stress_fraction(int dim,const double G[3][3],
        double parallel,double perpendicular,double *compression){
    if(dim<1||dim>3||!(parallel>=0)||!(perpendicular>=0)||!isfinite(parallel)||!isfinite(perpendicular))return -1;
    double trace=0,norm=0;
    for(int a=0;a<dim;a++)for(int b=0;b<dim;b++){
        if(!isfinite(G[a][b]))return -1;norm+=G[a][b]*G[a][b];if(a==b)trace+=G[a][b];
    }
    double negative=fmax(0,-trace);*compression=norm>0?fmin(1,negative*negative/norm):0;
    double denominator=parallel+perpendicular;
    return denominator>0?(perpendicular+*compression*parallel)/denominator:0;
}
#undef GFS_STRESS_HD
#endif
