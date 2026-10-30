#include "../gfs_kepler_vacuum_stencil.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

int main(void)
{
    GfsKeplerVacuumStencil s = {0};
    size_t degree1 = 0, degree2 = 0;
    assert(gfs_kepler_vacuum_stencil_build(
        256, 4.0, 0.5, 2.0, 1, 0.0, 0.0, &s));
    assert(s.computational_n == 258);
    assert(s.n_fluid == 48240);
    assert(s.n_sites == 908);
    assert(s.n_faces == 1280);
    assert(s.n_outside_target_box == 88);
    assert(fabs(s.dx - 1.0 / 64.0) < 1e-15);

    for (size_t k = 0; k < s.n_sites; ++k) {
        const GfsKeplerVacuumSite *site = &s.sites[k];
        const GfsKeplerVacuumIncidence *a = &s.incidences[4 * k];
        assert(site->stable_id == site->grid_j * s.computational_n +
                                  site->grid_i);
        assert(fabs(site->x - (-2.0 +
                    ((double)(site->grid_i - 1) + 0.5) * s.dx)) < 1e-14);
        assert(fabs(site->y - (-2.0 +
                    ((double)(site->grid_j - 1) + 0.5) * s.dx)) < 1e-14);
        if (site->incidence_count == 1) {
            ++degree1;
        } else {
            const GfsKeplerVacuumIncidence *b = a + 1;
            assert(site->incidence_count == 2);
            assert((int)a->nx * (int)b->nx + (int)a->ny * (int)b->ny == 0);
            ++degree2;
        }
        for (size_t q = 0; q < site->incidence_count; ++q) {
            const GfsKeplerVacuumIncidence *f = a + q;
            int fi = f->fluid_id % s.computational_n;
            int fj = f->fluid_id / s.computational_n;
            assert((f->nx == -1 || f->nx == 1) ? f->ny == 0 :
                   ((f->ny == -1 || f->ny == 1) && f->nx == 0));
            assert(site->grid_i == fi + f->nx);
            assert(site->grid_j == fj + f->ny);
        }
    }
    assert(degree1 == 536);
    assert(degree2 == 372);

    gfs_kepler_vacuum_stencil_destroy(&s);
    assert(s.sites == NULL && s.incidences == NULL && s.n_sites == 0);
    assert(!gfs_kepler_vacuum_stencil_build(
        256, 4.0, 2.0, 0.5, 1, 0.0, 0.0, &s));

    /* The original box puts the disk edge on the periodic seam. Its 88
     * first-layer vacuum sites outside the fundamental box could overlap
     * periodic material images, so the periodic constructor must reject it. */
    assert(!gfs_kepler_vacuum_stencil_build_periodic(
        256, 4.0, 0.5, 2.0, 1, 0.0, 0.0, &s));

    /* Expand the box while preserving dx and therefore the annulus sample. */
    assert(gfs_kepler_vacuum_stencil_build_periodic(
        320, 5.0, 0.5, 2.0, 1, 0.0, 0.0, &s));
    assert(s.n_fluid == 48240);
    assert(s.n_sites == 908);
    assert(s.n_faces == 1280);
    assert(s.n_outside_target_box == 0);
    assert(fabs(s.dx - 1.0 / 64.0) < 1e-15);
    gfs_kepler_vacuum_stencil_destroy(&s);
    return 0;
}
