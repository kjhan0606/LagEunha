#ifndef GFS_MATERIAL_STAR_ENTROPY_H
#define GFS_MATERIAL_STAR_ENTROPY_H
#include "gfs_material_pressure_jacobian.h"
#ifdef __CUDACC__
#define GFS_STAR_ENTROPY_HD static __host__ __device__ inline
#else
#define GFS_STAR_ENTROPY_HD static inline
#endif
/* Opt-in material-star candidate, common to dimensions 1--3. Nonnegative
 * face entropy production is sufficient, not necessary for an accurate
 * high-order scheme. This limiter is not yet a qualified hydro default.
 * Force, volume and TOTAL-energy exchange must use the returned pair as
 * -A p n, A u, -A p u. No independent heat correction is permitted. */
GFS_STAR_ENTROPY_HD double gfs_material_star_entropy_pair(double pl,double pr,
 double rl,double rr,double vl,double vr,double p,double u){
 return (rl/pl)*(pl-p)*(u-vl)+(rr/pr)*(pr-p)*(vr-u);
}
GFS_STAR_ENTROPY_HD int gfs_material_star_entropy_limit(int dim,
 GfsFvPrimitive l,GfsFvPrimitive r,const double normal[3],double gamma,
 double*p,double*u,double*theta){
 double vl=0,vr=0;for(int d=0;d<dim;d++){vl+=normal[d]*l.u[d];vr+=normal[d]*r.u[d];}
 double high=gfs_material_star_entropy_pair(l.p,r.p,l.rho,r.rho,vl,vr,*p,*u);
 *theta=1;if(high>=0)return 1;
 GfsMaterialFaceResponse low;if(gfs_material_face_response(dim,l,r,normal,gamma,&low)!=1)return 0;
 double dp=*p-low.pstar,du=*u-low.ustar,tl=l.rho/l.p,tr=r.rho/r.p;
 double C=gfs_material_star_entropy_pair(l.p,r.p,l.rho,r.rho,vl,vr,low.pstar,low.ustar);
 double scale=tl*fabs((l.p-low.pstar)*(low.ustar-vl))+tr*fabs((r.p-low.pstar)*(vr-low.ustar));
 if(C< -128*DBL_EPSILON*scale)return 0;
 if(C<0)C=0; /* Roundoff of the first-order certificate, not a heat floor. */
 double A=(tr-tl)*dp*du;
 double B=tl*((l.p-low.pstar)*du-dp*(low.ustar-vl))+tr*(-(r.p-low.pstar)*du-dp*(vr-low.ustar));
 double q=0;
 if(C==0){if(B>0&&A<0)q=-B/A;}
 else if(A==0){if(!(B<0))return 0;q=-C/B;}
 else{double D=B*B-4*A*C;if(!(D>=0))return 0;double s=sqrt(D);q=B>=0?(-B-s)/(2*A):2*C/(-B+s);}
 if(!isfinite(q)||q<0||q>1+128*DBL_EPSILON)return 0;q=fmin(q,1.0);
 *p=low.pstar+q*dp;*u=low.ustar+q*du;*theta=q;
 return isfinite(*p)&&isfinite(*u)&&*p>=0;
}
/* Differentiate with physical owner states fixed while reconstructed side
 * pressures vary. Weight correction must use this derivative when enabled.
 * Other response fields retain raw Riemann values and MUST NOT be used as
 * the impedance/Jacobian of the limited exchange. */
GFS_STAR_ENTROPY_HD int gfs_material_star_entropy_jacobian(int dim,
 GfsFvPrimitive owner_l,GfsFvPrimitive owner_r,
 GfsFvPrimitive side_l,GfsFvPrimitive side_r,const double normal[3],double gamma,
 GfsMaterialFaceResponse*star,GfsMaterialPressureJacobian*jac,double*theta){
 if(gfs_material_face_response(dim,side_l,side_r,normal,gamma,star)!=1)return 0;
 GfsMaterialPressureJacobian raw;
 if(!gfs_material_pressure_jacobian(side_l.rho,side_l.p,side_r.rho,side_r.p,star->pstar,gamma,&raw))return 0;
 double ph=star->pstar,uh=star->ustar;
 if(!gfs_material_star_entropy_limit(dim,owner_l,owner_r,normal,gamma,&star->pstar,&star->ustar,theta))return 0;
 *jac=raw;if(*theta==1)return 1;
 if(*theta==0){jac->dp_left=jac->dp_right=jac->du_left=jac->du_right=0;return 1;}
 GfsMaterialFaceResponse low;
 if(gfs_material_face_response(dim,owner_l,owner_r,normal,gamma,&low)!=1)return 0;
 double vl=0,vr=0;for(int d=0;d<dim;d++){vl+=normal[d]*owner_l.u[d];vr+=normal[d]*owner_r.u[d];}
 double dp=ph-low.pstar,du=uh-low.ustar,tl=owner_l.rho/owner_l.p,tr=owner_r.rho/owner_r.p;
 double sp=-tl*(star->ustar-vl)-tr*(vr-star->ustar);
 double su=tl*(owner_l.p-star->pstar)-tr*(owner_r.p-star->pstar);
 double den=sp*dp+su*du;
 if(den==0||!isfinite(den))return 0; /* Multiple root: no false derivative. */
 double dtl=-(*theta)*(sp*raw.dp_left+su*raw.du_left)/den;
 double dtr=-(*theta)*(sp*raw.dp_right+su*raw.du_right)/den;
 jac->dp_left=(*theta)*raw.dp_left+dp*dtl;jac->du_left=(*theta)*raw.du_left+du*dtl;
 jac->dp_right=(*theta)*raw.dp_right+dp*dtr;jac->du_right=(*theta)*raw.du_right+du*dtr;
 return isfinite(jac->dp_left)&&isfinite(jac->dp_right)&&isfinite(jac->du_left)&&isfinite(jac->du_right);
}
#undef GFS_STAR_ENTROPY_HD
#endif
