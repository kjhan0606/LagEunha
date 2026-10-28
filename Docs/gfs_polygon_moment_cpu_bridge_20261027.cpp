#include "gfs_polygon_moment_20261027.h"

extern "C" int gfs_polygon_moment_batch(
    int n,const int *offset,const double *x1,const double *y1,const double *x2,
    const double *y2,const double *mass,double *area,double *cx,double *cy,
    double *inertia,int *status)
{
    if(n<=0||!offset||!x1||!y1||!x2||!y2||!mass||!area||!cx||!cy||!inertia||!status)
        return -1;
    for(int i=0;i<n;++i) {
        gfs_polygon_moment p=gfs_polygon_moments(x1,y1,x2,y2,offset[i],offset[i+1],mass[i]);
        area[i]=p.area;cx[i]=p.cx;cy[i]=p.cy;inertia[i]=p.inertia;status[i]=p.status;
    }
    return 0;
}
