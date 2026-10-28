#ifndef EUNHA_AFFINE_FACE_TRANSPORT_H
#define EUNHA_AFFINE_FACE_TRANSPORT_H

/* Experimental 2D straight-face transport, not wired into production RHS.
 * z in [-A/2,A/2], normal mass flux j(z)=j0+j1*z,
 * position c+t*z, velocity v+g*z. Input: A,j0,j1,cx,cy,tx,ty,vx,vy,gx,gy.
 * Output: mass, first spatial moments x/y, momentum x/y, kinetic energy.
 * Caller supplies signed mass-flux convention and compatible thermal/pressure
 * terms. No positivity, Riemann solver, geometry or weight evolution is implied.
 */
#ifdef __CUDACC__
#define AFT_HD __host__ __device__
#else
#define AFT_HD
#endif
static AFT_HD inline void affine_face_transport(const double *a, double *f) {
    const double length=a[0], j0=a[1], j1=a[2];
    const double m=length*j0;
    const double second=length*length*length/12.0;
    f[0]=m;
    f[1]=m*a[3]+second*j1*a[5];
    f[2]=m*a[4]+second*j1*a[6];
    f[3]=m*a[7]+second*j1*a[9];
    f[4]=m*a[8]+second*j1*a[10];
    f[5]=0.5*m*(a[7]*a[7]+a[8]*a[8])
        +second*(j1*(a[7]*a[9]+a[8]*a[10])
        +0.5*j0*(a[9]*a[9]+a[10]*a[10]));
}
#undef AFT_HD
#endif
