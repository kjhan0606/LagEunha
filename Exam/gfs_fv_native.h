#ifndef GFS_FV_NATIVE_H
#define GFS_FV_NATIVE_H
#include <stdint.h>
#include "gfs_fv_ale.h"
enum {GFS_FV_PERIODIC=0,GFS_FV_REFLECT=1,GFS_FV_OUTFLOW=2};
typedef struct {
    double q[3],W,Q[5],mu,m0,W0,residual_integral;
    uint64_t id;
} GfsFvCell;
typedef struct {double x[3],weight;} GfsFvPoint;
typedef struct {
    int j,mirror,point,count;
    double r[3],normal[3],area,centroid[3],distance;
} GfsFvNativeFace;
typedef struct {
    double V,centroid[3]; /* centroid relative to own generator */
    GfsFvNativeFace *face;GfsFvPoint *point;
    int faces,face_capacity,points,point_capacity;
} GfsFvGeometry;
typedef struct {
    int dim,n,bc[3];double lo[3],L[3],gamma;
    GfsFvCell *cell;GfsFvGeometry *geometry;
    double geometry_seconds,maximum_volume_closure;
    int maximum_candidates;
} GfsFvMesh;
typedef struct {
    double dQ[5],qdot[3],z,Vdot,residual;
} GfsFvRate;
typedef struct {
    int newton,cg,backtracks,hll_points,points;
    double wave_dt,admissible_dt,root_residual,mass_link_defect;
    double boundary_flux[5],reconstruction_seconds,fan_seconds,solve_seconds,flux_seconds;
    double kinetic_rate,kinetic_pressure_work,stress_thermal_rate,viscous_dt;
    double entropy_rate,entropy_boundary_flux,entropy_dissipation,entropy_balance_defect;
} GfsFvRhsStats;
typedef struct GfsFvWorkspace GfsFvWorkspace;
int gfs_fv_geometry(GfsFvMesh *mesh);
void gfs_fv_geometry_free(GfsFvMesh *mesh);
GfsFvWorkspace *gfs_fv_workspace_create(void);
void gfs_fv_workspace_free(GfsFvWorkspace *workspace);
int gfs_fv_native_rhs(GfsFvWorkspace*,GfsFvMesh*,GfsFvRate*,GfsFvRhsStats*);
int gfs_fv_native_velocity(GfsFvWorkspace*,GfsFvMesh*,double (*velocity)[3]);
#endif
