#include <math.h>
#include <stdio.h>
#include "kepler_affine2d_stable.h"

int main(void)
{
    const double x[] = {3., 3.0001, 3.0001, 3.};
    const double y[] = {3., 3., 3.0001, 3.0001};
    const double width = x[1] - x[0];
    const double expected_area = width * width;
    const double expected_moment = expected_area * expected_area / 12.;
    KeplerPolygonMoments2D p;
    if (!kepler_polygon_moments2d(x, y, 4, &p)) return 1;
    double area_error = fabs(p.area / expected_area - 1.);
    double xx_error = fabs(p.sxx / expected_moment - 1.);
    double yy_error = fabs(p.syy / expected_moment - 1.);
    printf("relative area %.9g, sxx %.9g, syy %.9g\n", area_error, xx_error, yy_error);
    return area_error < 1e-10 && xx_error < 1e-10 && yy_error < 1e-10 ? 0 : 1;
}
