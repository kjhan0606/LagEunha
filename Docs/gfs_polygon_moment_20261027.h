#ifndef GFS_POLYGON_MOMENT_20261027_H
#define GFS_POLYGON_MOMENT_20261027_H

#include <math.h>

#if defined(__CUDACC__)
#define GFS_POLY_HD __host__ __device__
#else
#define GFS_POLY_HD
#endif

/* Endpoints are consistently oriented polygon edges in coordinates relative
 * to the mesh generator. The returned centroid is an offset from that
 * generator; inertia assumes uniform density and is about the centroid.
 * The second pass centers coordinates before the polar moment sum. */
struct gfs_polygon_moment {
    double area,cx,cy,polar_area_moment,inertia,closure,perimeter;
    int status;
};

GFS_POLY_HD static inline gfs_polygon_moment gfs_polygon_moments(
    const double *x1,const double *y1,const double *x2,const double *y2,
    int begin,int end,double mass)
{
    gfs_polygon_moment out={0.,0.,0.,0.,0.,0.,0.,-1};
    if (!x1 || !y1 || !x2 || !y2 || begin<0 || end-begin<3 ||
        !isfinite(mass) || !(mass>0.)) return out;
    double cross_sum=0.,mx=0.,my=0.,sx=0.,sy=0.,perimeter=0.;
    for(int k=begin;k<end;++k) {
        double ax=x1[k],ay=y1[k],bx=x2[k],by=y2[k];
        if(!isfinite(ax)||!isfinite(ay)||!isfinite(bx)||!isfinite(by)) return out;
        double dx=bx-ax,dy=by-ay,len=sqrt(dx*dx+dy*dy);
        if(!isfinite(len)||!(len>0.)) return out;
        double cross=ax*by-bx*ay;
        cross_sum+=cross;mx+=(ax+bx)*cross;my+=(ay+by)*cross;
        sx+=dx;sy+=dy;perimeter+=len;
    }
    out.perimeter=perimeter;
    out.closure=sqrt(sx*sx+sy*sy);
    if(!isfinite(out.closure)||out.closure>1e-10*perimeter) { out.status=-2;return out; }
    double signed_area=.5*cross_sum;
    if(!isfinite(signed_area)||signed_area==0.) { out.status=-3;return out; }
    out.cx=mx/(6.*signed_area);out.cy=my/(6.*signed_area);
    if(!isfinite(out.cx)||!isfinite(out.cy)) { out.status=-3;return out; }
    double moment=0.;
    for(int k=begin;k<end;++k) {
        double ax=x1[k]-out.cx,ay=y1[k]-out.cy;
        double bx=x2[k]-out.cx,by=y2[k]-out.cy;
        double cross=ax*by-bx*ay;
        moment+=cross*(ax*ax+ax*bx+bx*bx+ay*ay+ay*by+by*by);
    }
    out.area=fabs(signed_area);
    out.polar_area_moment=(signed_area>0. ? moment : -moment)/12.;
    out.inertia=mass*out.polar_area_moment/out.area;
    if(!isfinite(out.inertia)||!(out.inertia>0.)) { out.status=-4;return out; }
    out.status=0;
    return out;
}

#undef GFS_POLY_HD
#endif
