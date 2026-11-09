#ifndef GFS_FIXED_CYLINDER_CONTACT_H
#define GFS_FIXED_CYLINDER_CONTACT_H
#include "gfs_material_fixed_solid2d.h"
#include "gfs_material_fixed_cylinder_data.h"
#include <stdint.h>
/* Same polygon as geometry, no circular surrogate. Positive signed distance
 * means exterior; the nearest supporting plane is selected only inside S. */
static inline int gfs_fixed_cylinder_outside_box(double x,double y){return x<2.||x>3.||y<1.5||y>2.5;}
static inline double gfs_fixed_cylinder_gap(double x,double y,int*face,double*n){
 double best=-HUGE_VAL;
 for(int k=0;k<GFS_FIXED_CYLINDER_N;k++){
  int j=(k+1)%GFS_FIXED_CYLINDER_N;
  double ex=gfs_fixed_cylinder_x[j]-gfs_fixed_cylinder_x[k],ey=gfs_fixed_cylinder_y[j]-gfs_fixed_cylinder_y[k],A=hypot(ex,ey);
  double d=-gfs_solid_cross(ex,ey,x-gfs_fixed_cylinder_x[k],y-gfs_fixed_cylinder_y[k])/A;
  if(d>best){best=d;*face=k;n[0]=ey/A;n[1]=-ex/A;}
 }
 return best;
}
static inline int gfs_fixed_cylinder_contact(double*x,double*y,double*vx,double*vy){
 if(!isfinite(*x)||!isfinite(*y)||!isfinite(*vx)||!isfinite(*vy))return 0;
 if(gfs_fixed_cylinder_outside_box(*x,*y))return 1;
 int face;double n[2],gap=gfs_fixed_cylinder_gap(*x,*y,&face,n);
 if(gap>=0)return 1;
 *x=fma(-2*gap,n[0],*x);*y=fma(-2*gap,n[1],*y);
 double vn=*vx*n[0]+*vy*n[1];
 if(vn<0){*vx-=vn*n[0];*vy-=vn*n[1];}
 /* Only the representable exterior condition, never a geometric epsilon. */
 double qn[2];int q;
 if(gfs_fixed_cylinder_gap(*x,*y,&q,qn)<0){
  if(fabs(n[0])>=fabs(n[1]))*x=nextafter(*x,n[0]>0?HUGE_VAL:-HUGE_VAL);
  else *y=nextafter(*y,n[1]>0?HUGE_VAL:-HUGE_VAL);
 }
 return gfs_fixed_cylinder_gap(*x,*y,&q,qn)>=0;
}
static inline uint64_t gfs_fixed_cylinder_signature(void){
 uint64_t h=UINT64_C(14695981039346656037);
 for(int k=0;k<GFS_FIXED_CYLINDER_N;k++){
  double p[2]={gfs_fixed_cylinder_x[k],gfs_fixed_cylinder_y[k]};
  const unsigned char*b=(const unsigned char*)p;
  for(unsigned j=0;j<sizeof(p);j++){h^=b[j];h*=UINT64_C(1099511628211);}
 }
 return h;
}
#endif
