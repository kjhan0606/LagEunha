#ifndef GFS_CIRCLE_NATIVE_FACES_H
#define GFS_CIRCLE_NATIVE_FACES_H
#include "gfs_circle_geometry.h"
#ifdef __CUDACC__
#define GFS_CNF_HD __host__ __device__
#else
#define GFS_CNF_HD
#endif

/* A curved-wall record carries endpoint-integrated normals, not arc length
 * as a straight-face area. The native W0 wall branch uses only this normal
 * integral for work; circular swept-volume absolute values use the arc law.
 */
#define GFS_CIRCLE_ARC_SUPPORT (-2000000000)

GFS_CNF_HD static inline int gfs_circle_native_metrics(
    const double *x,const double *y,int n,double cx,double cy,double radius,
    GfsCircleMetrics *out){
    if(!x || !y || !out || n<3 || !(radius>0.) || !isfinite(radius))return 0;
    double cross=0.,angle=0.,nx=0.,ny=0.;
    for(int k=0;k<n;k++){
        int j=(k+1)%n;
        GfsCirclePoint a={x[k]+cx,y[k]+cy},b={x[j]+cx,y[j]+cy};
        if(!isfinite(a.x)||!isfinite(a.y)||!isfinite(b.x)||!isfinite(b.y))return 0;
        if(a.x==b.x && a.y==b.y)continue;
        GfsCircleSegment s;
        if(gfs_circle_clip_segment(a,b,radius,&s)){
            cross+=gfs_circle_cross(s.a,s.b);
            nx+=s.b.y-s.a.y;ny-=s.b.x-s.a.x;
            if(s.u0>0.)angle+=gfs_circle_angle(a,s.a);
            if(s.u1<1.)angle+=gfs_circle_angle(s.b,b);
        }else angle+=gfs_circle_angle(a,b);
    }
    out->area=.5*(cross+radius*radius*angle);
    out->wall_arc_length=radius*angle;out->wall_nx=-nx;out->wall_ny=-ny;
    return isfinite(out->area)&&isfinite(out->wall_arc_length)&&isfinite(nx)&&isfinite(ny);
}

GFS_CNF_HD static inline int gfs_circle_native_emit(
    double x,double y,int support,double *ox,double *oy,int *os,int *count,int cap){
    int n=*count;
    if(n && ox[n-1]==x && oy[n-1]==y){os[n-1]=support;return 1;}
    if(n>=cap)return 0;
    ox[n]=x;oy[n]=y;os[n]=support;*count=n+1;return 1;
}

/* Reuse the caller's existing polygon/scratch buffers. Original artificial
 * box supports (-1..-4) are retained if they still intersect the disk, so the
 * native adaptive search must expand instead of accepting a truncated cell.
 * Circumscribed tangent supports (<=-1000) must have zero resolved length.
 * Return -1 for invalid/capacity/unresolved topology, otherwise new count.
 * changed=0 means the original polygon is wholly inside; no copy occurred.
 */
GFS_CNF_HD static inline int gfs_circle_native_outline(
    const double *x,const double *y,const int *support,int n,double cx,double cy,
    double radius,double *ox,double *oy,int *os,int cap,int *changed,double *area){
    if(!x || !y || !support || !ox || !oy || !os || !changed || !area ||
       n<3 || n>cap || cap<4 || !(radius>0.) || !isfinite(radius))return -1;
    int inside=1;
    for(int k=0;k<n;k++){
        double a=x[k]+cx,b=y[k]+cy;
        if(!isfinite(a)||!isfinite(b))return -1;
        if(a*a+b*b>radius*radius)inside=0;
        if(support[k]<=-1000)inside=0;
    }
    *changed=0;*area=-1.;
    if(inside)return n;
    GfsCircleMetrics metrics;
    if(!gfs_circle_native_metrics(x,y,n,cx,cy,radius,&metrics) || !(metrics.area>0.))return -1;
    double close=32.*DBL_EPSILON*fmax(radius,fabs(cx)+fabs(cy));
    int m=0;
    for(int k=0;k<n;k++){
        int j=(k+1)%n;
        GfsCirclePoint a={x[k]+cx,y[k]+cy},b={x[j]+cx,y[j]+cy};
        GfsCircleSegment s;
        if(!gfs_circle_clip_segment(a,b,radius,&s))continue;
        if(support[k]<=-1000)return -1;
        /* Interpolate in owner-relative coordinates. Preserve uncut endpoint
         * bits and support-plane accuracy instead of adding/subtracting the
         * much larger global circle-center offset at every real vertex.
         */
        GfsCirclePoint local_a={x[k],y[k]},local_b={x[j],y[j]};
        GfsCirclePoint clipped_a=gfs_circle_lerp(local_a,local_b,s.u0);
        GfsCirclePoint clipped_b=gfs_circle_lerp(local_a,local_b,s.u1);
        double ax=clipped_a.x,ay=clipped_a.y,bx=clipped_b.x,by=clipped_b.y;
        if(m && hypot(ox[m-1]-ax,oy[m-1]-ay)<=close)os[m-1]=support[k];
        else if(!gfs_circle_native_emit(ax,ay,support[k],ox,oy,os,&m,cap))return -1;
        if(!gfs_circle_native_emit(bx,by,GFS_CIRCLE_ARC_SUPPORT,ox,oy,os,&m,cap))return -1;
    }
    if(!m){
        if(fabs(metrics.area-acos(-1.)*radius*radius)>128.*DBL_EPSILON*radius*radius)return -1;
        ox[0]=radius-cx;oy[0]=-cy;
        ox[1]=-cx;oy[1]=radius-cy;
        ox[2]=-radius-cx;oy[2]=-cy;
        ox[3]=-cx;oy[3]=-radius-cy;
        for(int k=0;k<4;k++)os[k]=GFS_CIRCLE_ARC_SUPPORT;
        m=4;
    }else if(m>1 && hypot(ox[m-1]-ox[0],oy[m-1]-oy[0])<=close)m--;
    if(m<2)return -1;
    *changed=1;*area=metrics.area;
    return m;
}

GFS_CNF_HD static inline double gfs_circle_native_arc_angle(
    GfsCirclePoint a,GfsCirclePoint b){
    double angle=gfs_circle_angle(a,b);
    if(angle<0.)angle+=2.*acos(-1.);
    return angle;
}

/* Exact |v.n| integral along a CCW circle arc, splitting at at most two
 * zeros of the normal speed. This prevents within-arc sign cancellation
 * from underestimating the geometric CFL sweep.
 */
GFS_CNF_HD static inline double gfs_circle_native_arc_sweep(
    GfsCirclePoint a,GfsCirclePoint b,double radius,double vx,double vy){
    if(vx==0. && vy==0.)return 0.;
    double pi=acos(-1.),lo=atan2(a.y,a.x),hi=lo+gfs_circle_native_arc_angle(a,b);
    double root=atan2(vy,vx)+.5*pi;
    double previous=vx*a.y-vy*a.x,total=0.;
    for(int k=-4;k<=4;k++){
        double t=root+k*pi;
        if(t>lo && t<hi){
            double value=radius*(vx*sin(t)-vy*cos(t));
            total+=fabs(value-previous);previous=value;
        }
    }
    total+=fabs(vx*b.y-vy*b.x-previous);
    double net=fabs(vx*(b.y-a.y)-vy*(b.x-a.x));
    return fmax(total,net);
}

#undef GFS_CNF_HD
#endif
