#ifndef GFS_CIRCLE_GEOMETRY_H
#define GFS_CIRCLE_GEOMETRY_H

/* Analytic intersection of one convex CCW polygon with a fixed circle.
 * All coordinates are relative to the circle center. Intended as a shared
 * host/CUDA building block, not yet connected to the production solver.
 * Straight material interfaces retain their original neighbor/support IDs.
 * Circular wall normal integral follows exactly from contour closure.
 */
#include <math.h>
#include <float.h>
#ifdef __CUDACC__
#define GFS_CIRCLE_HD __host__ __device__
#else
#define GFS_CIRCLE_HD
#endif

typedef struct { double x,y; } GfsCirclePoint;
typedef struct {
    double u0,u1;
    GfsCirclePoint a,b;
} GfsCircleSegment;
typedef struct {
    double area, wall_arc_length, wall_nx, wall_ny;
} GfsCircleMetrics;

GFS_CIRCLE_HD static inline double gfs_circle_cross(GfsCirclePoint a,GfsCirclePoint b){
    return a.x*b.y-a.y*b.x;
}
GFS_CIRCLE_HD static inline double gfs_circle_angle(GfsCirclePoint a,GfsCirclePoint b){
    return atan2(gfs_circle_cross(a,b),a.x*b.x+a.y*b.y);
}
GFS_CIRCLE_HD static inline GfsCirclePoint gfs_circle_lerp(
    GfsCirclePoint a,GfsCirclePoint b,double u){
    if(u==0.)return a;
    if(u==1.)return b;
    GfsCirclePoint p={a.x+u*(b.x-a.x),a.y+u*(b.y-a.y)};
    return p;
}

/* Clip a straight interface to the disk. A tangent point has zero length. */
GFS_CIRCLE_HD static inline int gfs_circle_clip_segment(
    GfsCirclePoint a,GfsCirclePoint b,double radius,GfsCircleSegment *out){
    if(!out || !isfinite(radius) || !isfinite(a.x) || !isfinite(a.y) ||
       !isfinite(b.x) || !isfinite(b.y))return 0;
    double dx=b.x-a.x,dy=b.y-a.y;
    double A=dx*dx+dy*dy;
    if(!(A>0.) || !(radius>0.))return 0;
    double B=a.x*dx+a.y*dy;
    double C=a.x*a.x+a.y*a.y-radius*radius;
    double D=fma(-A,C,B*B);
    /* Bound evaluation roundoff in A, B, C and the discriminant. An
     * unresolved tangency has zero segment measure; do not manufacture
     * a sqrt(epsilon)-length face from a roundoff-positive discriminant.
     * The bound scales with the actual arithmetic, not a physical length.
     */
    double error=16.*DBL_EPSILON*(B*B+fabs(A*C)+
                    A*(a.x*a.x+a.y*a.y+radius*radius));
    if(!(D>error))return 0;
    double root=sqrt(D);
    double q=-B-copysign(root,B);
    if(q==0.)return 0;
    double v0=q/A,v1=C/q;
    if(v0>v1){double t=v0;v0=v1;v1=t;}
    double u0=fmax(0.,v0),u1=fmin(1.,v1);
    if(!(u1>u0))return 0;
    out->u0=u0;out->u1=u1;
    out->a=gfs_circle_lerp(a,b,u0);
    out->b=gfs_circle_lerp(a,b,u1);
    return 1;
}

/* Sector integration replaces only the portions of the original contour
 * outside the circle. Internal cells use the shoelace fast path, with no
 * square roots or inverse trigonometric functions.
 */
GFS_CIRCLE_HD static inline int gfs_circle_polygon_metrics(
    const GfsCirclePoint *vertices,int n,double radius,GfsCircleMetrics *out){
    if(!vertices || !out || n<3 || !(radius>0.) || !isfinite(radius))return 0;
    int inside=1;
    double radius2=radius*radius;
    for(int k=0;k<n;k++){
        if(!isfinite(vertices[k].x) || !isfinite(vertices[k].y))return 0;
        if(vertices[k].x*vertices[k].x+vertices[k].y*vertices[k].y>radius2)inside=0;
    }
    double twice_area=0.,angle=0.,line_nx=0.,line_ny=0.;
    for(int k=0;k<n;k++){
        GfsCirclePoint a=vertices[k],b=vertices[(k+1)%n];
        if(inside){twice_area+=gfs_circle_cross(a,b);continue;}
        if(a.x==b.x && a.y==b.y)continue;
        GfsCircleSegment segment;
        if(gfs_circle_clip_segment(a,b,radius,&segment)){
            twice_area+=gfs_circle_cross(segment.a,segment.b);
            line_nx+=segment.b.y-segment.a.y;
            line_ny-=segment.b.x-segment.a.x;
            if(segment.u0>0.)angle+=gfs_circle_angle(a,segment.a);
            if(segment.u1<1.)angle+=gfs_circle_angle(segment.b,b);
        }else angle+=gfs_circle_angle(a,b);
    }
    out->area=.5*(twice_area+radius2*angle);
    out->wall_arc_length=radius*angle;
    out->wall_nx=-line_nx;out->wall_ny=-line_ny;
    return isfinite(out->area) && isfinite(out->wall_arc_length) &&
           isfinite(out->wall_nx) && isfinite(out->wall_ny);
}

#undef GFS_CIRCLE_HD
#endif
