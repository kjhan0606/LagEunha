#ifndef GFS_VISCOSITY_IMPLICIT_H
#define GFS_VISCOSITY_IMPLICIT_H
/* Positive central viscosity on a canonical undirected graph. Coefficients
 * are lagged; endpoint compression signs are solved. BE is first order.
 * x and v use interleaved dim*n doubles. Walls have j=-1 and outward n.
 * Internal n points from i to j and must be parallel to x_j-x_i.
 */
#ifdef __cplusplus
extern "C" {
#endif
typedef struct {
    int i,j;
    double normal[3],coefficient;
} GfsViscosityEdge;
typedef struct GfsViscosityWorkspace GfsViscosityWorkspace;
typedef struct {
    int iterations,active_iterations,active_edges;
    double relative_equation_residual,absolute_equation_residual;
    double energy_defect,heat,min_heat,max_heat_over_U;
    double momentum_defect[3],angular_defect[3];
    double wall_impulse[3],wall_torque[3];
    double solve_seconds;
} GfsViscosityStats;
GfsViscosityWorkspace *gfs_viscosity_workspace_create(void);
void gfs_viscosity_workspace_free(GfsViscosityWorkspace *w);
/* Returns 0 on success. A failure leaves v and U unchanged. */
int gfs_viscosity_be(GfsViscosityWorkspace *w,int dim,int n,int ne,
    const GfsViscosityEdge *edges,const double *mass,const double *x,
    double *v,double *U,double interval,GfsViscosityStats *stats);
#ifdef __cplusplus
}
#endif
#endif
