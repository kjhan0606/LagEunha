#include <stddef.h>
#include <omp.h>
#include "moment_quadrupole_gravity.h"
int mq_threads(int n) {
    int actual=0;omp_set_dynamic(0);omp_set_num_threads(n);
    #pragma omp parallel reduction(max:actual)
    actual=omp_get_num_threads();
    return actual;
}
int mq_batch(size_t n,const double *in,double *out) {
    int invalid=0;
    #pragma omp parallel for schedule(static) reduction(|:invalid)
    for(size_t k=0;k<n;k++) invalid |= mq_source(in+18*k,out+12*k);
    return invalid;
}
