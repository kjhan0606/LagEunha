#ifndef GFS_MATERIAL_RADIAL_IMAGE_H
#define GFS_MATERIAL_RADIAL_IMAGE_H
#include <math.h>
#include <float.h>
#include "gfs_power_predicate2d.h"
#ifdef __CUDACC__
#define GFS_RADIAL_HD static __host__ __device__ inline
#else
#define GFS_RADIAL_HD static inline
#endif
/* Radial reflection through a fixed sphere (circle in 2D).
 * Geometry: g(x)=center+(2R-r)n. Fluid free-slip reflection is S=I-2nn^T.
 * Coordinate derivative is J=(2R/r-1)(I-nn^T)-nn^T, not S.
 * J is symmetric: this same map pulls a ghost force back to its owner.
 * Neither map declares the resulting polygonal power boundary stationary.
 */
typedef struct {
    double point[3], normal[3], tangent_scale;
} GfsMaterialRadialImage;
GFS_RADIAL_HD int gfs_material_radial_image(int dim,const double *center,
        double radius,const double *x,GfsMaterialRadialImage *image){
    if(dim<1||dim>3||!center||!x||!image||!(radius>0)||!isfinite(radius))return 0;
    double d[3]={0,0,0},r=0;
    for(int k=0;k<dim;k++){
        if(!isfinite(center[k])||!isfinite(x[k]))return 0;
        d[k]=x[k]-center[k];r=hypot(r,d[k]);
    }
    if(!(r>0)||!isfinite(r))return 0;
    double reflected=2*radius-r,scale=2*(radius/r)-1;
    if(!isfinite(reflected)||!isfinite(scale))return 0;
    for(int k=0;k<3;k++){
        image->normal[k]=k<dim?d[k]/r:0;
        image->point[k]=k<dim?center[k]+image->normal[k]*reflected:0;
        if(!isfinite(image->point[k]))return 0;
    }
    image->tangent_scale=scale;return 1;
}
GFS_RADIAL_HD int gfs_material_radial_vector(int dim,const GfsMaterialRadialImage *image,
        const double *v,int coordinate_derivative,double *out){
    if(dim<1||dim>3||!image||!v||!out)return 0;
    double vn=0,norm=0,scale=coordinate_derivative?image->tangent_scale:1;
    if(!isfinite(scale))return 0;
    for(int k=0;k<dim;k++){
        if(!isfinite(v[k])||!isfinite(image->normal[k]))return 0;
        vn+=v[k]*image->normal[k];norm+=image->normal[k]*image->normal[k];
    }
    if(!isfinite(vn)||fabs(norm-1)>64*DBL_EPSILON)return 0;
    for(int k=0;k<dim;k++){
        /* Explicit fluid branch preserves the existing reflection arithmetic. */
        out[k]=coordinate_derivative?scale*v[k]-(scale+1)*vn*image->normal[k]
                                    :v[k]-2*vn*image->normal[k];
        if(!isfinite(out[k]))return 0;
    }
    return 1;
}
/* Same reflected point, stored as owner plus a local displacement. Retaining
 * the low part avoids losing the image normal when the owner approaches R. */
GFS_RADIAL_HD int gfs_material_radial_split_point(int dim,const double*center,
 double radius,const double*x,double*high,double*low){
 GfsMaterialRadialImage a;
 if(!high||!low||!gfs_material_radial_image(dim,center,radius,x,&a))return 0;
 double r=0;for(int k=0;k<dim;k++)r=hypot(r,x[k]-center[k]);
 double displacement=2*(radius-r);
 for(int k=0;k<dim;k++){
  GfsPowerDD p=gfs_pdd_sum(x[k],displacement*a.normal[k]);
  high[k]=p.hi;low[k]=p.lo;
  if(!isfinite(p.hi)||!isfinite(p.lo))return 0;
 }
 return 1;
}
/* Normal speed of the owner's own equal-weight image face. Differentiate
 * n.(X-center)=R before dividing: no inverse wall gap or speed cap remains.
 * Sign is outward from the fluid cell toward its reflected site. */
GFS_RADIAL_HD int gfs_material_radial_self_face_speed(int dim,const double*center,
 double radius,const double*x,const double*relative_centroid,const double*v,double*speed){
 GfsMaterialRadialImage a;
 if(!relative_centroid||!v||!speed||!gfs_material_radial_image(dim,center,radius,x,&a))return 0;
 double r=0,cn=0,vn=0;
 for(int k=0;k<dim;k++){
  if(!isfinite(relative_centroid[k])||!isfinite(v[k]))return 0;
  r=hypot(r,x[k]-center[k]);cn+=relative_centroid[k]*a.normal[k];vn+=v[k]*a.normal[k];
 }
 if(r==radius)return 0;
 double dot=0;
 for(int k=0;k<dim;k++)dot+=(relative_centroid[k]-cn*a.normal[k])*(v[k]-vn*a.normal[k]);
 *speed=(r>radius?dot:-dot)/r;
 return isfinite(*speed);
}
#undef GFS_RADIAL_HD
#endif
