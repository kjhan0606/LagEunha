#ifndef GFS_MATERIAL_MAP12_H
#define GFS_MATERIAL_MAP12_H
#include <math.h>
#ifdef __CUDACC__
#define GFS_M12_HD __host__ __device__
#else
#define GFS_M12_HD
#endif
#define GFS_MATERIAL12_N 12
/* Compatible finite material kinematics and symmetric stress-work algebra.
 * x and v share the SAME fixed reference basis. Caller may center x and add
 * that fixed origin to physical position afterwards. No time integrator,
 * gravity model, thermal reset, viscosity or remap is supplied here.
 */
typedef struct {int n;double x[2][12],v[2][12];} GfsMaterialMap12;
typedef struct {double position[2],velocity[2],F[4],J,grad_v[4];} GfsMaterialPoint12;
GFS_M12_HD static inline int gfs_material12_point(const GfsMaterialMap12 *m,
    const double *N,const double *Nx,const double *Ny,GfsMaterialPoint12 *p){
    if(m->n<1 || m->n>12)return 0;
    for(int k=0;k<2;k++)p->position[k]=p->velocity[k]=0;
    for(int k=0;k<4;k++)p->F[k]=p->grad_v[k]=0;
    double V[4]={0};
    for(int a=0;a<m->n;a++)for(int k=0;k<2;k++){
        p->position[k]+=m->x[k][a]*N[a];p->velocity[k]+=m->v[k][a]*N[a];
        p->F[2*k]+=m->x[k][a]*Nx[a];p->F[2*k+1]+=m->x[k][a]*Ny[a];
        V[2*k]+=m->v[k][a]*Nx[a];V[2*k+1]+=m->v[k][a]*Ny[a];
    }
    p->J=p->F[0]*p->F[3]-p->F[1]*p->F[2];if(!(p->J>0) || !isfinite(p->J))return 0;
    double inv[4]={p->F[3]/p->J,-p->F[1]/p->J,-p->F[2]/p->J,p->F[0]/p->J};
    for(int k=0;k<2;k++){p->grad_v[2*k]=V[2*k]*inv[0]+V[2*k+1]*inv[2];p->grad_v[2*k+1]=V[2*k]*inv[1]+V[2*k+1]*inv[3];}return 1;
}
/* Use an already evaluated point when pressure/body force depends on its
 * physical position/J. This avoids a duplicate geometry evaluation on GPU.
 */
GFS_M12_HD static inline int gfs_material12_stress_at_point(int n,
    const GfsMaterialPoint12 *p,const double *N,const double *Nx,const double *Ny,
    double rho0,double weight,const double sigma[3],double gx,double gy,
    double force[2][12],double *Udot,double *gravity_power){
    if(n<1 || n>12 || !(rho0>0) || !(weight>=0) || !(p->J>0) || !isfinite(p->J))return 0;
    double inv[4]={p->F[3]/p->J,-p->F[1]/p->J,-p->F[2]/p->J,p->F[0]/p->J},dm=weight*rho0,dV=weight*p->J;
    for(int a=0;a<n;a++){
        double dx=Nx[a]*inv[0]+Ny[a]*inv[2],dy=Nx[a]*inv[1]+Ny[a]*inv[3];
        force[0][a]+=-dV*(sigma[0]*dx+sigma[1]*dy)+dm*N[a]*gx;
        force[1][a]+=-dV*(sigma[1]*dx+sigma[2]*dy)+dm*N[a]*gy;
    }
    *Udot+=dV*(sigma[0]*p->grad_v[0]+sigma[1]*(p->grad_v[1]+p->grad_v[2])+sigma[2]*p->grad_v[3]);
    *gravity_power+=dm*(p->velocity[0]*gx+p->velocity[1]*gy);return 1;
}
GFS_M12_HD static inline int gfs_material12_stress_point(const GfsMaterialMap12 *m,
    const double *N,const double *Nx,const double *Ny,double rho0,double weight,
    const double sigma[3],double gx,double gy,double force[2][12],double *Udot,double *gravity_power){
    GfsMaterialPoint12 p;if(!gfs_material12_point(m,N,Nx,Ny,&p))return 0;
    return gfs_material12_stress_at_point(m->n,&p,N,Nx,Ny,rho0,weight,sigma,gx,gy,force,Udot,gravity_power);
}
/* Reference traction = (J sigma F^-T)n0. ratio=dS/dS0 for unit n0.
 * Tangent/normal gradients may differ across cells; shared position/velocity
 * traces and a single hybrid traction still determine pair work/torque.
 */
GFS_M12_HD static inline int gfs_material12_piola_traction(const GfsMaterialPoint12 *p,
    const double sigma[3],double nx,double ny,double traction[2],double *ratio){
    if(!(p->J>0))return 0;double ax=p->F[3]*nx-p->F[2]*ny,ay=-p->F[1]*nx+p->F[0]*ny;
    traction[0]=sigma[0]*ax+sigma[1]*ay;traction[1]=sigma[1]*ax+sigma[2]*ay;*ratio=hypot(ax,ay);return *ratio>0 && isfinite(*ratio);
}
GFS_M12_HD static inline int gfs_material12_mass_point(int n,const double *N,
    double dm,double M[12][12]){
    if(n<1 || n>12 || !(dm>=0))return 0;for(int a=0;a<n;a++)for(int b=0;b<n;b++)M[a][b]+=dm*N[a]*N[b];return 1;
}
GFS_M12_HD static inline double gfs_material12_kinetic(const GfsMaterialMap12 *m,const double M[12][12]){
    double K=0;for(int i=0;i<m->n;i++)for(int j=0;j<m->n;j++)K+=.5*M[i][j]*(m->v[0][i]*m->v[0][j]+m->v[1][i]*m->v[1][j]);return K;
}
typedef struct {
    double mass,center[2],momentum[2],kinetic,internal_kinetic;
    double intrinsic_angular,polar_inertia,rigid_rotation_kinetic,resolved_remainder;
} GfsMaterialMoments12;
/* constant[a] represents the function1; e.g. ten ones and two zeros for P3.
 * Use centered quadratic forms for internal K/I/L, avoiding bulk subtraction.
 * These are actual represented-field moments, never an independent reservoir.
 */
GFS_M12_HD static inline int gfs_material12_moments(const GfsMaterialMap12 *m,
    const double M[12][12],const double constant[12],GfsMaterialMoments12 *out){
    if(m->n<1 || m->n>12)return 0;out->mass=0;double q[2]={0},p[2]={0};
    for(int i=0;i<m->n;i++)for(int j=0;j<m->n;j++){
        double r=constant[i]*M[i][j];out->mass+=r*constant[j];for(int k=0;k<2;k++){q[k]+=r*m->x[k][j];p[k]+=r*m->v[k][j];}
    }
    if(!(out->mass>0))return 0;double xc[2][12],vc[2][12];
    for(int k=0;k<2;k++){out->center[k]=q[k]/out->mass;out->momentum[k]=p[k];for(int i=0;i<m->n;i++){xc[k][i]=m->x[k][i]-out->center[k]*constant[i];vc[k][i]=m->v[k][i]-p[k]/out->mass*constant[i];}}
    out->internal_kinetic=out->intrinsic_angular=out->polar_inertia=0;
    for(int i=0;i<m->n;i++)for(int j=0;j<m->n;j++){
        out->internal_kinetic+=.5*M[i][j]*(vc[0][i]*vc[0][j]+vc[1][i]*vc[1][j]);out->polar_inertia+=M[i][j]*(xc[0][i]*xc[0][j]+xc[1][i]*xc[1][j]);out->intrinsic_angular+=M[i][j]*(xc[0][i]*vc[1][j]-xc[1][i]*vc[0][j]);
    }
    if(!(out->polar_inertia>0))return 0;out->kinetic=gfs_material12_kinetic(m,M);out->rigid_rotation_kinetic=.5*out->intrinsic_angular*out->intrinsic_angular/out->polar_inertia;out->resolved_remainder=out->internal_kinetic-out->rigid_rotation_kinetic;
    return isfinite(out->resolved_remainder);
}
#undef GFS_M12_HD
#endif
