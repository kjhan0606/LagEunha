#ifndef GFS_KEPLER_VACUUM_STENCIL_H
#define GFS_KEPLER_VACUUM_STENCIL_H

#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct {
    int grid_i, grid_j;       /* index in the guarded computational lattice */
    int stable_id;            /* row-major ID in that full lattice */
    double x, y;              /* initial generator position */
    size_t incidence_count;
} GfsKeplerVacuumSite;

typedef struct {
    int fluid_id;             /* row-major ID in the guarded full lattice */
    signed char nx, ny;       /* initial cardinal normal, fluid -> generator */
} GfsKeplerVacuumIncidence;

typedef struct {
    int target_n, guard, computational_n;
    double dx;
    size_t n_fluid, n_sites, n_faces, n_outside_target_box;
    GfsKeplerVacuumSite *sites;
    /* Four fixed slots per site keep this small setup allocation simple. */
    GfsKeplerVacuumIncidence *incidences;
} GfsKeplerVacuumStencil;

static inline void gfs_kepler_vacuum_stencil_destroy(
    GfsKeplerVacuumStencil *stencil)
{
    if (!stencil) return;
    free(stencil->sites);
    free(stencil->incidences);
    *stencil = (GfsKeplerVacuumStencil){0};
}

/* Build the W=0, one-cardinal-layer geometry for a cell-centered annulus.
 * `box` is the width of the unguarded n x n target lattice. The guard extends
 * the computational lattice without changing target coordinates or spacing.
 * This returns geometry metadata only; it does not create solver particles
 * or prescribe later-time generator motion. `out` must be zero-initialized or
 * destroyed with gfs_kepler_vacuum_stencil_destroy before reuse. */
static inline int gfs_kepler_vacuum_stencil_build(
    int n, double box, double r_in, double r_out, int guard,
    double cx, double cy, GfsKeplerVacuumStencil *out)
{
    static const int dir_x[4] = {-1, 1, 0, 0};
    static const int dir_y[4] = {0, 0, -1, 1};
    GfsKeplerVacuumStencil result = {0};
    unsigned char *fluid = NULL;
    int *site_at = NULL;
    size_t *written = NULL;
    size_t n_grid, n_fluid = 0, n_sites = 0, n_faces = 0;
    int n_comp, i, j, d;

    if (!out || n <= 0 || guard < 1 || guard > (INT_MAX - n) / 2 ||
        !(box > 0.0) || !(r_in >= 0.0) || !(r_out > r_in) ||
        !isfinite(box) || !isfinite(r_in) || !isfinite(r_out) ||
        !isfinite(cx) || !isfinite(cy))
        return 0;

    n_comp = n + 2 * guard;
    if ((size_t)n_comp > SIZE_MAX / (size_t)n_comp) return 0;
    n_grid = (size_t)n_comp * (size_t)n_comp;
    if (n_grid > (size_t)INT_MAX || (size_t)n > SIZE_MAX / (size_t)n)
        return 0;

    result.target_n = n;
    result.guard = guard;
    result.computational_n = n_comp;
    result.dx = box / (double)n;
    fluid = (unsigned char *)calloc(n_grid, sizeof(*fluid));
    site_at = (int *)malloc(n_grid * sizeof(*site_at));
    if (!fluid || !site_at) goto fail;
    for (size_t k = 0; k < n_grid; ++k) site_at[k] = -1;

    for (j = 0; j < n; ++j) {
        double y = cy - 0.5 * box + ((double)j + 0.5) * result.dx;
        for (i = 0; i < n; ++i) {
            double x = cx - 0.5 * box + ((double)i + 0.5) * result.dx;
            double r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
            if (r2 >= r_in * r_in && r2 <= r_out * r_out) {
                size_t gi = (size_t)(i + guard);
                size_t gj = (size_t)(j + guard);
                fluid[gj * (size_t)n_comp + gi] = 1;
                ++n_fluid;
            }
        }
    }
    if (n_fluid == 0 || n_fluid > SIZE_MAX / 4) goto fail;

    /* First pass: unique missing cardinal coordinates and face counts. */
    result.sites = (GfsKeplerVacuumSite *)calloc(4 * n_fluid,
                                                  sizeof(*result.sites));
    if (!result.sites) goto fail;
    for (j = 0; j < n; ++j) {
        for (i = 0; i < n; ++i) {
            size_t gi = (size_t)(i + guard);
            size_t gj = (size_t)(j + guard);
            size_t fluid_id = gj * (size_t)n_comp + gi;
            if (!fluid[fluid_id]) continue;
            for (d = 0; d < 4; ++d) {
                int qx = i + dir_x[d], qy = j + dir_y[d];
                int gx = qx + guard, gy = qy + guard;
                size_t qid;
                int site_index;
                if (gx < 0 || gx >= n_comp || gy < 0 || gy >= n_comp)
                    goto fail;
                qid = (size_t)gy * (size_t)n_comp + (size_t)gx;
                if (fluid[qid]) continue;
                site_index = site_at[qid];
                if (site_index < 0) {
                    GfsKeplerVacuumSite *site;
                    if (n_sites >= 4 * n_fluid || n_sites > (size_t)INT_MAX)
                        goto fail;
                    site_index = (int)n_sites++;
                    site_at[qid] = site_index;
                    site = &result.sites[site_index];
                    site->grid_i = gx;
                    site->grid_j = gy;
                    site->stable_id = (int)qid;
                    site->x = cx - 0.5 * box + ((double)qx + 0.5) * result.dx;
                    site->y = cy - 0.5 * box + ((double)qy + 0.5) * result.dx;
                    if (qx < 0 || qx >= n || qy < 0 || qy >= n)
                        ++result.n_outside_target_box;
                }
                if (result.sites[site_index].incidence_count >= 4) goto fail;
                ++result.sites[site_index].incidence_count;
                ++n_faces;
            }
        }
    }
    if (n_sites == 0 || n_sites > SIZE_MAX / 4) goto fail;
    result.incidences = (GfsKeplerVacuumIncidence *)calloc(
        4 * n_sites, sizeof(*result.incidences));
    written = (size_t *)calloc(n_sites, sizeof(*written));
    if (!result.incidences || !written) goto fail;

    /* Second pass: attach stable fluid IDs and outward lattice normals. */
    for (j = 0; j < n; ++j) {
        for (i = 0; i < n; ++i) {
            size_t gi = (size_t)(i + guard);
            size_t gj = (size_t)(j + guard);
            size_t fluid_id = gj * (size_t)n_comp + gi;
            if (!fluid[fluid_id]) continue;
            for (d = 0; d < 4; ++d) {
                int gx = i + dir_x[d] + guard;
                int gy = j + dir_y[d] + guard;
                size_t qid = (size_t)gy * (size_t)n_comp + (size_t)gx;
                int site_index = site_at[qid];
                size_t slot;
                GfsKeplerVacuumSite *site;
                GfsKeplerVacuumIncidence *incidence;
                if (fluid[qid]) continue;
                if (site_index < 0 || (size_t)site_index >= n_sites) goto fail;
                site = &result.sites[site_index];
                slot = written[site_index]++;
                if (slot >= site->incidence_count || slot >= 4) goto fail;
                slot += 4 * (size_t)site_index;
                incidence = &result.incidences[slot];
                incidence->fluid_id = (int)fluid_id;
                incidence->nx = (signed char)dir_x[d];
                incidence->ny = (signed char)dir_y[d];
            }
        }
    }

    /* The audited target has one face per site or two orthogonal faces. */
    for (size_t s = 0; s < n_sites; ++s) {
        GfsKeplerVacuumSite *site = &result.sites[s];
        GfsKeplerVacuumIncidence *a = &result.incidences[4 * s];
        if (site->incidence_count == 0 || site->incidence_count > 2) goto fail;
        if (site->incidence_count == 2) {
            GfsKeplerVacuumIncidence *b = a + 1;
            if ((int)a->nx * (int)b->nx + (int)a->ny * (int)b->ny != 0)
                goto fail;
        }
    }

    result.n_fluid = n_fluid;
    result.n_sites = n_sites;
    result.n_faces = n_faces;
    free(fluid);
    free(site_at);
    free(written);
    *out = result;
    return 1;

fail:
    free(fluid);
    free(site_at);
    free(written);
    gfs_kepler_vacuum_stencil_destroy(&result);
    return 0;
}

/* Build an annular stencil for a square periodic domain centered at (cx,cy).
 * Require a full cell-width of vacuum between the outer fluid support and
 * the periodic seam. The geometric site audit is retained as a second guard:
 * every one-layer vacuum generator must remain inside the fundamental box,
 * where periodic copies of material cells cannot overlap it. */
static inline int gfs_kepler_vacuum_stencil_build_periodic(
    int n, double box, double r_in, double r_out, int guard,
    double cx, double cy, GfsKeplerVacuumStencil *out)
{
    GfsKeplerVacuumStencil result = {0};
    double dx;

    if (!out || n <= 0 || !(box > 0.0) || !isfinite(box) ||
        !isfinite(r_out) || !(r_out >= 0.0))
        return 0;
    dx = box / (double)n;
    if (!(r_out < 0.5 * box - dx))
        return 0;
    if (!gfs_kepler_vacuum_stencil_build(
            n, box, r_in, r_out, guard, cx, cy, &result))
        return 0;
    if (result.n_outside_target_box != 0) {
        gfs_kepler_vacuum_stencil_destroy(&result);
        return 0;
    }
    *out = result;
    return 1;
}

#endif
