#ifndef GPU_VALIDATION_CHECK_H
#define GPU_VALIDATION_CHECK_H
#include <math.h>

/* Same mixed tolerance convention as the original ax/energy checks.
 * Nonfinite values must never pass through unordered floating comparisons. */
static inline int gpu_validation_mismatch(double actual, double reference,
    double relative_tolerance, double absolute_tolerance, double *relative_error)
{
    if (!isfinite(actual) || !isfinite(reference)) {
        *relative_error = INFINITY;
        return 1;
    }
    double difference = fabs(actual-reference);
    double scale = fmax(fabs(actual),fabs(reference)) + 1.0e-30;
    *relative_error = difference/scale;
    return difference > absolute_tolerance && *relative_error > relative_tolerance;
}
#endif
