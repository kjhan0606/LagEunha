#ifndef GFS_MATERIAL_FREE_BUDGET2D_H
#define GFS_MATERIAL_FREE_BUDGET2D_H
#include <math.h>
#include <stddef.h>
/* Frozen input to exactly the existing finite-predictor inequality.
 * This is a certificate for reusing an UNPROJECTED RHS, not a new limiter. */
typedef struct {double dx,dy,dvx,dvy,dw,dq;} GfsMaterialFreePair2D;
typedef struct {double dx,dy;} GfsMaterialFreeSupport2D;
static inline int gfs_material_free_pairs_support_valid2d(int count,const GfsMaterialFreePair2D*p,const GfsMaterialFreeSupport2D*support,double h,double beta){
 if(count<1||!p||!(h>0)||!isfinite(h)||!(beta>0&&beta<1))return 0;
 for(int k=0;k<count;k++){
  double dx=p[k].dx+h*p[k].dvx,dy=p[k].dy+h*p[k].dvy;
  double dw=p[k].dw+h*p[k].dq,d2=dx*dx+dy*dy;
  if(support){double sx=support[k].dx,sy=support[k].dy;d2=2*(sx*dx+sy*dy)-sx*sx-sy*sy;}
  if(!(d2>0)||!isfinite(d2)||!isfinite(dw)||fabs(dw)>nextafter(beta*d2,0.0))return 0;
 }
 return 1;
}
static inline int gfs_material_free_pairs_valid2d(int count,const GfsMaterialFreePair2D*p,double h,double beta){
 return gfs_material_free_pairs_support_valid2d(count,p,NULL,h,beta);
}
#endif
