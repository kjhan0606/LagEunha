#ifndef GFS_VISCOSITY_IMPLICIT_H
#define GFS_VISCOSITY_IMPLICIT_H
/* Positive central viscosity on a canonical undirected graph. Coefficients
 * are lagged; endpoint compression signs are solved. BE is first order.
 * x and v use interleaved dim*n doubles. Walls have j=-1 and outward n.
 * Internal n points from i to the selected image of j and must be parallel
 * to x_j + image_shift - x_i. Zero image_shift gives the ordinary domain.
 */
#ifdef __cplusplus
extern "C" {
#endif
typedef struct {
    int i,j;
    double normal[3],coefficient;
    double image_shift[3]; /* fixed periodic translation, no velocity shift */
} GfsViscosityEdge;
typedef struct GfsViscosityWorkspace GfsViscosityWorkspace;
typedef struct {
    int iterations,active_iterations,active_edges;
    double relative_equation_residual,absolute_equation_residual;
    double energy_defect,heat,min_heat,max_heat_over_U;
    double momentum_defect[3],angular_defect[3];
    double wall_impulse[3],wall_torque[3];
    double periodic_torque[3]; /* angular-momentum flux across periodic cuts */
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
