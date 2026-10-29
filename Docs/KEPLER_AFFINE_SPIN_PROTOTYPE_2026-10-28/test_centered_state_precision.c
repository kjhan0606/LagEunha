#include <math.h>
#include <stdio.h>
#include "kepler_affine2d_stable.h"

int main(void)
{
    const double width = 1e-4, m = 1., c = 3., omega = 3.;
    const double moment = m * width * width / 12.;
    KeplerAffineConserved2D raw = {0};
    KeplerAffineStateLagrangian2D centered = {0};
    KeplerAffineDerived2D a, b;
    raw.mass = m; raw.hx = raw.hy = m*c;
    raw.qxx = raw.qyy = m*c*c + moment; raw.qxy = m*c*c;
    raw.rxy = omega*moment; raw.ryx = -omega*moment;
    raw.etotal = 1.;
    centered.cx = centered.cy = c;
    centered.sxx = centered.syy = moment;
    centered.bxy = omega*moment; centered.byx = -omega*moment;
    centered.etotal = 1.;
    if (!kepler_affine_derive2d(&raw, &a) || !kepler_affine_state_derive2d(&centered, m, 0., 0., &b)) return 1;
    const double raw_error = fabs(a.gyx/omega-1.);
    const double centered_error = fabs(b.gyx/omega-1.);
    printf("raw gradient relative error %.9g, centered relative error %.9g\n", raw_error, centered_error);
    return raw_error > 1e-8 && centered_error < 1e-14 ? 0 : 1;
}
