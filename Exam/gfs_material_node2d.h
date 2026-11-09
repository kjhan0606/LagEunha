#ifndef GFS_MATERIAL_NODE2D_H
#define GFS_MATERIAL_NODE2D_H

#ifdef __CUDACC__
#define GFS_NODE_HD static __host__ __device__ inline
#else
#define GFS_NODE_HD static inline
#endif
#include "gfs_exact_riemann_ale.h"
/* Half-edge contributions at a shared material vertex. Normals point OUT
 * of each incident cell. Length is half the cell edge length. No mass flux.
 * The single vertex velocity balances all exact piston-pressure tractions.
 * This local operator does not provide ownership, mesh persistence or PLM.
 */
typedef struct {double n[2],length,rho,p,u[2];} GfsMaterialHalfEdge;
GFS_NODE_HD int gfs_material_node_system(int count,const GfsMaterialHalfEdge *e,
 double gamma,const double v[2],double F[2],double M[3],double *scale){
 F[0]=F[1]=M[0]=M[1]=M[2]=*scale=0;
 for(int k=0;k<count;k++){
  double s=(e[k].u[0]-v[0])*e[k].n[0]+(e[k].u[1]-v[1])*e[k].n[1],p,d;
  if(!(e[k].length>0)||!gfs_exact_piston(e[k].rho,e[k].p,gamma,s,&p,&d))return 0;
  double nx=e[k].n[0],ny=e[k].n[1],l=e[k].length;
  if(fabs(nx*nx+ny*ny-1)>1e-10)return 0;
  F[0]+=l*p*nx;F[1]+=l*p*ny;
  M[0]+=l*d*nx*nx;M[1]+=l*d*nx*ny;M[2]+=l*d*ny*ny;
  *scale+=l*(e[k].p+p);
 }
 return isfinite(F[0])&&isfinite(F[1])&&isfinite(*scale);
}
GFS_NODE_HD int gfs_material_node_solve(int count,const GfsMaterialHalfEdge *e,
 double gamma,int free_x,int free_y,double v[2],int *iterations,double *relative_residual){
 if(count<1)return 0;
 for(int it=0;it<64;it++){
  double F[2],M[3],scale;
  if(!gfs_material_node_system(count,e,gamma,v,F,M,&scale))return 0;
  double r=hypot(free_x?F[0]:0,free_y?F[1]:0);
  *iterations=it;*relative_residual=r/fmax(scale,1e-300);
  if(r<=2e-13*fmax(scale,1e-300))return 1;
  double d[2]={0,0};
  if(free_x&&free_y){
   double det=M[0]*M[2]-M[1]*M[1];
   if(!(det>0)||!isfinite(det))return 0;
   d[0]=(M[2]*F[0]-M[1]*F[1])/det;
   d[1]=(M[0]*F[1]-M[1]*F[0])/det;
  }else if(free_x){if(!(M[0]>0))return 0;d[0]=F[0]/M[0];}
  else if(free_y){if(!(M[2]>0))return 0;d[1]=F[1]/M[2];}
  else return 1;
  int accepted=0;
  for(int back=0;back<32;back++){
   double step=ldexp(1.,-back),trial[2]={v[0]+step*d[0],v[1]+step*d[1]};
   double f[2],m[3],sc;
   if(gfs_material_node_system(count,e,gamma,trial,f,m,&sc)&&hypot(free_x?f[0]:0,free_y?f[1]:0)<r){
    v[0]=trial[0];v[1]=trial[1];accepted=1;break;
   }
  }
  if(!accepted)return 0;
 }
 return 0;
}
/* Cell force and energy gain use exactly the SAME solved vertex velocity.
 * heat_excess is the irreversible part relative to -P dV for a cell-mean
 * velocity, and is nonnegative for the monotone exact piston curve.
 */
GFS_NODE_HD int gfs_material_halfedge_rate(const GfsMaterialHalfEdge *e,
 double gamma,const double v[2],double force[2],double *energy,double *heat_excess){
 double slip=(v[0]-e->u[0])*e->n[0]+(v[1]-e->u[1])*e->n[1],p,d;
 if(!gfs_exact_piston(e->rho,e->p,gamma,-slip,&p,&d))return 0;
 force[0]=-e->length*p*e->n[0];force[1]=-e->length*p*e->n[1];
 *energy=force[0]*v[0]+force[1]*v[1];
 *heat_excess=-e->length*(p-e->p)*slip;
 return isfinite(*energy)&&isfinite(*heat_excess);
}
#undef GFS_NODE_HD
#endif
