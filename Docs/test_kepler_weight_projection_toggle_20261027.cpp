#include "gfs_pressure_neutral_projection_20261027.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

int main() {
    const unsigned char active[3] = {1, 1, 0};
    const double gradient[3] = {2., 1., 99.};
    const double original[3] = {3., 4., 5.};
    double rate[3], alpha = -1.;
    std::memcpy(rate, original, sizeof rate);
    if (setenv("GFS_GPU_VOLW2_PROJECTION", "0", 1)) return 1;
    int status = gfs_pressure_neutral_projection(rate, gradient, active, 3, &alpha);
    if (status || alpha != 0. || std::memcmp(rate, original, sizeof rate)) return 2;
    std::memcpy(rate, original, sizeof rate);
    if (setenv("GFS_GPU_VOLW2_PROJECTION", "1", 1)) return 3;
    if (gfs_pressure_neutral_projection(rate, gradient, active, 3, &alpha)
        || !std::isfinite(alpha) || std::fabs(gradient[0]*rate[0]+gradient[1]*rate[1]) > 1e-14) return 4;
    if (unsetenv("GFS_GPU_VOLW2_PROJECTION")) return 5;
    std::memcpy(rate, original, sizeof rate);
    if (gfs_pressure_neutral_projection(rate, gradient, active, 3, &alpha)
        || !std::isfinite(alpha) || std::fabs(gradient[0]*rate[0]+gradient[1]*rate[1]) > 1e-14) return 6;
    if (setenv("GFS_GPU_VOLW2_PROJECTION", "typo", 1)) return 7;
    if (gfs_pressure_neutral_projection(rate, gradient, active, 3, &alpha) != -5) return 8;
    std::puts("PASS: raw-rate switch, default projection, and invalid mode");
    return 0;
}
