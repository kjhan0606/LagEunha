#include <stddef.h>
#include <omp.h>
#include "weight_remap_point_energy.h"
int wp_threads(int n) {
    int actual=0; omp_set_dynamic(0); omp_set_num_threads(n);
    #pragma omp parallel reduction(max:actual)
    actual=omp_get_num_threads();
    return actual;
}
int wp_batch(size_t n,const double *in,double *out) {
    #pragma omp parallel for schedule(static)
    for(size_t k=0;k<n;k++) wp_pair(in+9*k,out+6*k);
    return 0;
}
