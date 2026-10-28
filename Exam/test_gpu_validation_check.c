#include <assert.h>
#include <stdio.h>
#include "gpu_validation_check.h"

int main(void)
{
    double error;
    assert(!gpu_validation_mismatch(1.,1.,1e-6,1e-11,&error));
    assert(!gpu_validation_mismatch(0.,1e-12,1e-6,1e-11,&error));
    assert(!gpu_validation_mismatch(1.,1.+1e-7,1e-6,1e-11,&error));
    assert(gpu_validation_mismatch(1.,1.01,1e-6,1e-11,&error));
    assert(gpu_validation_mismatch(-1.,1.,1e-6,1e-11,&error));
    assert(gpu_validation_mismatch(NAN,1.,1e-6,1e-11,&error));
    assert(isinf(error));
    assert(gpu_validation_mismatch(1.,NAN,1e-6,1e-11,&error));
    assert(gpu_validation_mismatch(INFINITY,INFINITY,1e-6,1e-11,&error));
    assert(gpu_validation_mismatch(-INFINITY,0.,1e-6,1e-11,&error));
    /* Timestep comparison has no absolute floor that could hide a large
     * relative error in a small but positive timestep. */
    assert(gpu_validation_mismatch(1e-12,2e-12,1e-6,0.,&error));
    assert(!gpu_validation_mismatch(1e30,1e30,1e-6,0.,&error));
    assert(!gpu_validation_mismatch(0.,1e-9,1e-6,1e-8,&error));
    puts("GPU validation scalar checks: 12 cases passed");
    return 0;
}
