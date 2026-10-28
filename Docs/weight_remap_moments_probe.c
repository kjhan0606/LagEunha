#include <stddef.h>
#include <omp.h>
#include "weight_remap_moments.h"

int wm_set_threads(int n) {
    int actual=0;
    omp_set_dynamic(0); omp_set_num_threads(n);
    #pragma omp parallel reduction(max:actual)
    actual=omp_get_num_threads();
    return actual;
}
int wm_batch(size_t n, const double *in, double *out) {
    #pragma omp parallel for schedule(static)
    for (size_t k=0;k<n;k++) wm_atom(in+7*k,out+17*k);
    return 0;
}
int wm_lift_batch(size_t n, const double *in, double *out) {
    int invalid=0;
    #pragma omp parallel for schedule(static) reduction(|:invalid)
    for (size_t k=0;k<n;k++) invalid |= wm_lift_atom(in+16*k,out+17*k);
    return invalid;
}
