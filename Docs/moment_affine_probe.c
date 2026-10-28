#include <stddef.h>
#include <omp.h>
#include "moment_affine_recovery.h"
int ma_threads(int n) {
    int actual=0;omp_set_dynamic(0);omp_set_num_threads(n);
    #pragma omp parallel reduction(max:actual)
    actual=omp_get_num_threads();
    return actual;
}
int ma_batch(size_t n,const double *in,double *out) {
    int invalid=0;
    #pragma omp parallel for schedule(static) reduction(|:invalid)
    for(size_t k=0;k<n;k++) invalid |= ma_recover(in+17*k,out+21*k);
    return invalid;
}
