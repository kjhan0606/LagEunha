#ifndef GFS_ACOUSTIC_CFL_H
#define GFS_ACOUSTIC_CFL_H
#include "gfs_variational_core.h"

/* Rotation-invariant block-norm bound for the POSITIVE acoustic Hessian.
 * T_ij = sqrt(gamma P_i/V_i) |B_ij| / sqrt(m_j).
 * lambda_max(M^-1/2 B^T D B M^-1/2) <= max_row(T)*max_column(T).
 * Pressure-weighted geometric second derivatives are NOT bounded here.
 * Reciprocal interior faces are required for the incoming column sum. */
typedef struct {
    int dim;
    double root_D,inverse_root_mass,diagonal[3],row,column;
} GfsAcousticBound;
static inline int gfs_acoustic_bound_init(GfsAcousticBound *s,int dim,
        double mass,double volume,double U,double gamma)
{
    if(dim<1||dim>3||!(mass>0)||!(volume>0)||!(U>0)||!(gamma>1))return 0;
    s->dim=dim;s->root_D=sqrt(gamma*(gamma-1)*U)/volume;
    s->inverse_root_mass=1/sqrt(mass);s->row=s->column=0;
    for(int a=0;a<3;a++)s->diagonal[a]=0;
    return isfinite(s->root_D)&&isfinite(s->inverse_root_mass);
}
static inline int gfs_acoustic_bound_face(GfsAcousticBound *s,
        const GfsVariationalFace *face,double mass_j,double volume_j,
        double U_j,double gamma)
{
    if(face->dim!=s->dim||!(mass_j>0)||!(volume_j>0)||!(U_j>0))return 0;
    double bi2=0,bj2=0;
    for(int a=0;a<s->dim;a++){
        s->diagonal[a]-=face->bj[a];
        bi2+=face->bi[a]*face->bi[a];bj2+=face->bj[a]*face->bj[a];
    }
    s->row+=s->root_D*sqrt(bj2/mass_j);
    s->column+=sqrt(gamma*(gamma-1)*U_j*bi2)/volume_j*s->inverse_root_mass;
    return isfinite(s->row)&&isfinite(s->column);
}
static inline void gfs_acoustic_bound_finish(GfsAcousticBound *s)
{
    double norm=0;for(int a=0;a<s->dim;a++)norm+=s->diagonal[a]*s->diagonal[a];
    double self=s->root_D*sqrt(norm)*s->inverse_root_mass;
    s->row+=self;s->column+=self;
}
static inline double gfs_acoustic_bound_timestep(double courant,double row,double column)
{
    if(!(courant>0)||!(row>0)||!(column>0)||!isfinite(row)||!isfinite(column))return 0;
    return 2*courant/(sqrt(row)*sqrt(column));
}
#endif
