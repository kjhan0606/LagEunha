#ifndef LAGFORCE_WEIGHT_REMAP_MOMENTS_H
#define LAGFORCE_WEIGHT_REMAP_MOMENTS_H

/* Experimental positive-quadrature moment contribution for weight remapping.
 * Input: dm,x,y,vx,vy,e_thermal,phi_external at ONE spatial quadrature point.
 * Layout: m,qx,qy,Qxx,Qxy,Qyy,px,py,Dxx,Dxy,Dyx,Dyy,
 *         Pixx,Pixy,Piyy,E_hydro,E_potential.
 * D_ab = integral rho x_a v_b; Lz=Dxy-Dyx.
 * U=E_hydro-(Pixx+Piyy)/2: unresolved velocity variance is NOT heat.
 * These are raw volume moments, NOT the existing point-particle state.
 * Caller must provide a positive geometric partition and donor states.
 * No reconstruction, pressure force, gravity evolution, or restart conversion
 * is supplied here. Not connected to the production solver.
 */
#ifdef __CUDACC__
#define WM_INLINE static __host__ __device__ __forceinline__
#else
#define WM_INLINE static inline
#endif

WM_INLINE void wm_atom(const double *a, double *b) {
    double m=a[0], x=a[1], y=a[2], vx=a[3], vy=a[4];
    b[0]=m; b[1]=m*x; b[2]=m*y;
    b[3]=m*x*x; b[4]=m*x*y; b[5]=m*y*y;
    b[6]=m*vx; b[7]=m*vy;
    b[8]=m*x*vx; b[9]=m*x*vy; b[10]=m*y*vx; b[11]=m*y*vy;
    b[12]=m*vx*vx; b[13]=m*vx*vy; b[14]=m*vy*vy;
    b[15]=m*(a[5]+0.5*(vx*vx+vy*vy)); b[16]=m*a[6];
}

/* Evaluate the experimentally fitted positive density on a donor fragment.
 * Input: dV,x,y,Xi,Yi,vx,vy,rho_bar,thermal_energy_per_volume,
 *        phi(x),phi(Xi),lambda_x,lambda_y,lambda_L,lambda_phi,mean_lambda.
 * Lambda is already divided by the feature scale. mean_lambda is the
 * donor's mean of the scaled features dotted with the unscaled multiplier.
 * Returns nonzero for a nonpositive input state; caller must reject the
 * entire batch and must not consume its partial output.
 */
WM_INLINE int wm_lift_atom(const double *a, double *b) {
    double dx=a[1]-a[3],dy=a[2]-a[4];
    double feature=dx*a[11]+dy*a[12]+(dx*a[6]-dy*a[5])*a[13]
                  +(a[9]-a[10])*a[14]-a[15];
    double density=a[7]*(1.0+a[7]*feature);
    if (!(a[0]>0.0 && density>0.0 && a[8]>0.0)) return 1;
    double atom[7]={density*a[0],a[1],a[2],a[5],a[6],a[8]/density,a[9]};
    wm_atom(atom,b);
    return 0;
}
#undef WM_INLINE
#endif
