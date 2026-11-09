#ifndef GFS_MATERIAL_POTENTIAL2D_H
#define GFS_MATERIAL_POTENTIAL2D_H
#include "gfs_material_trace_csr2d.h"
/* Candidate geometry derivative only, not a gravity evolution operator.
 * Epot=sum m_i grad_phi.dot(volume_centroid_i), fixed m, uniform cell rho.
 * centroid_relative is relative to each generator. Directed faces must have
 * validated reverse partners; count each interior face once (i<j). Own-image
 * fixed walls have zero motion. Moving walls are not supported.
 * gx=dEpot/d(generator position), gW=dEpot/d(weight).
 */
static inline int gfs_material_potential_gradient2d(int n,int nf,
 const GfsTraceDirected2D *edge,const unsigned char *fixed_wall,
 const double *rho,const double *centroid_relative,const double grad_phi[2],
 double *gx,double *gW){
 if(n<1||nf<1||!edge||!rho||!centroid_relative||!grad_phi||!gx||!gW)return 0;
 if(!isfinite(grad_phi[0])||!isfinite(grad_phi[1]))return 0;
 for(int i=0;i<n;i++){
  if(!(rho[i]>0)||!isfinite(rho[i])||!isfinite(centroid_relative[2*i])||
     !isfinite(centroid_relative[2*i+1]))return 0;
  gx[2*i]=gx[2*i+1]=gW[i]=0;
 }
 for(int k=0;k<nf;k++){
  const GfsTraceDirected2D *e=edge+k;int i=e->i,j=e->j;
  int wall=fixed_wall?fixed_wall[k]:0;
  if(i<0||j<0||i>=n||j>=n||wall>1)return 0;
  if(wall){if(i!=j)return 0;continue;}
  if(i==j)return 0;if(i>j)continue;
  double d=hypot(e->dx,e->dy),A=e->area;
  if(!(d>0)||!isfinite(d)||!(A>0)||!isfinite(A)||!isfinite(e->cx)||!isfinite(e->cy))return 0;
  double nx=e->dx/d,ny=e->dy/d,tx=-ny,ty=nx;
  double pi=grad_phi[0]*(e->cx-centroid_relative[2*i])+grad_phi[1]*(e->cy-centroid_relative[2*i+1]);
  double pj=grad_phi[0]*(e->cx-e->dx-centroid_relative[2*j])+grad_phi[1]*(e->cy-e->dy-centroid_relative[2*j+1]);
  double psi=rho[i]*pi-rho[j]*pj;
  /* Integral of the product of affine potential and affine face speed.
   * A straight segment's central second moment is A^3/12. */
  double covariance=A*A*A*(rho[i]-rho[j])*(grad_phi[0]*tx+grad_phi[1]*ty)/(12*d);
  double ox=(e->cx-.5*e->dx)/d,oy=(e->cy-.5*e->dy)/d;
  gx[2*i]+=A*psi*(.5*nx+ox)+covariance*tx;
  gx[2*i+1]+=A*psi*(.5*ny+oy)+covariance*ty;
  gx[2*j]+=A*psi*(.5*nx-ox)-covariance*tx;
  gx[2*j+1]+=A*psi*(.5*ny-oy)-covariance*ty;
  gW[i]+=A*psi/(2*d);gW[j]-=A*psi/(2*d);
 }
 for(int i=0;i<n;i++)if(!isfinite(gx[2*i])||!isfinite(gx[2*i+1])||!isfinite(gW[i]))return 0;
 return 1;
}
#endif
