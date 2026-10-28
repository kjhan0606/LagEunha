/* Standalone OpenMP adapter for the experimental shared CPU/CUDA header. */
#include <stddef.h>
#include <omp.h>
#include "weight_remap_flux.h"

int wr_set_threads(int n) {
    int actual = 0;
    omp_set_dynamic(0);
    omp_set_num_threads(n);
    #pragma omp parallel reduction(max:actual)
    actual = omp_get_num_threads();
    return actual;
}

int wr_flux_batch(int mode, size_t n, const double *in, double *out) {
    if (mode != WR_DONOR && mode != WR_CENTERED && mode != WR_KINETIC_THERMAL_DONOR) return 1;
    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < n; ++k) {
        const double *a = in + 7*k;
        WrFlux f = wr_pair_flux(mode,a[0],a[1],a[2],a[3],a[4],a[5],a[6]);
        out[4*k] = f.mass; out[4*k+1] = f.px;
        out[4*k+2] = f.py; out[4*k+3] = f.energy;
    }
    return 0;
}
