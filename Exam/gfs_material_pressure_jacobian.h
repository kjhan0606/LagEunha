#ifndef GFS_MATERIAL_PRESSURE_JACOBIAN_H
#define GFS_MATERIAL_PRESSURE_JACOBIAN_H
#include "gfs_material_face_response.h"
#ifdef __CUDACC__
#define GFS_PRESSURE_JAC_HD static __host__ __device__ inline
#else
#define GFS_PRESSURE_JAC_HD static inline
#endif
/* Exact ideal-gas wave derivatives at a known non-vacuum contact pressure.
 * rho and normal velocity are held fixed. This is NOT a secant impedance.
 * Scalar normal Riemann problem: shared unchanged by dimensions 1, 2 and 3.
 * Diagnostic helper only; no production flux or pressure state is changed. */
GFS_PRESSURE_JAC_HD int gfs_material_wave_pressure_derivatives(double rho,
 double pressure,double star,double gamma,double *dstar,double *dside){
 if(!dstar||!dside||!(rho>0)||!(pressure>0)||!(star>0)||!(gamma>1)||
 !isfinite(rho)||!isfinite(pressure)||!isfinite(star)||!isfinite(gamma))return 0;
 double a,b;
 if(star>pressure){
  double beta=(gamma-1)/(gamma+1),den=star+beta*pressure;
  double root=sqrt(2/((gamma+1)*rho*den));
  a=root*(1-.5*(star-pressure)/den);
  b=-root*(1+.5*beta*(star-pressure)/den);
 }else{
  double c=sqrt(gamma*pressure/rho),power=pow(star/pressure,(gamma-1)/(2*gamma));
  a=c*power/(gamma*star);
  b=c/pressure*((power-1)/(gamma-1)-power/gamma);
 }
 if(!(a>0)||!isfinite(a)||!isfinite(b))return 0;
 *dstar=a;*dside=b;return 1;
}
typedef struct {double dp_left,dp_right,du_left,du_right;} GfsMaterialPressureJacobian;
GFS_PRESSURE_JAC_HD int gfs_material_pressure_jacobian(double rho_left,double p_left,
 double rho_right,double p_right,double pstar,double gamma,
 GfsMaterialPressureJacobian *out){
 double al,bl,ar,br;
 if(!out||!gfs_material_wave_pressure_derivatives(rho_left,p_left,pstar,gamma,&al,&bl)||
 !gfs_material_wave_pressure_derivatives(rho_right,p_right,pstar,gamma,&ar,&br))return 0;
 double sum=al+ar;GfsMaterialPressureJacobian q;
 if(!(sum>0)||!isfinite(sum))return 0;
 q.dp_left=-bl/sum;q.dp_right=-br/sum;
 q.du_left=ar*q.dp_left;q.du_right=-al*q.dp_right;
 if(!isfinite(q.dp_left)||!isfinite(q.dp_right)||!isfinite(q.du_left)||!isfinite(q.du_right))return 0;
 *out=q;return 1;
}
#undef GFS_PRESSURE_JAC_HD
#endif
