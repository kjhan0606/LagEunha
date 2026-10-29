#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include HEADER

int main(int argc, char **argv)
{
    if (argc != 2) return 2;
    double w = strtod(argv[1], NULL), x[] = {3., 3.+w, 3.+w, 3.};
    double y[] = {3., 3., 3.+w, 3.+w};
    double actual_w = x[1]-x[0], area = actual_w*actual_w, moment = area*area/12.;
    KeplerPolygonMoments2D p;
    if (!kepler_polygon_moments2d(x, y, 4, &p)) return 1;
    printf("width=%.17g area_relative=%.9g moment_relative=%.9g\n", actual_w,
           fabs(p.area/area-1.), fabs(p.sxx/moment-1.));
    return 0;
}
