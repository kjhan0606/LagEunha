#include <stddef.h>
#include <omp.h>
#include "moment_gravity_step.h"
int ms_threads(int n) {
    int actual=0;omp_set_dynamic(0);omp_set_num_threads(n);
    #pragma omp parallel reduction(max:actual)
    actual=omp_get_num_threads();
    return actual;
}
int ms_batch(size_t n,const double *in,double *out) {
    int invalid=0;
    #pragma omp parallel for schedule(static) reduction(|:invalid)
    for(size_t k=0;k<n;k++) invalid |= ms_step(in+19*k,out+17*k);
    return invalid;
}
