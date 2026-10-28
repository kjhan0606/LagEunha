#include "gfs_polygon_moments_reference_20261028.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>

static void require_close(const char *name, long double got,
		long double expected, long double tolerance)
{
	const long double error = std::abs(got - expected);
	if (!(error <= tolerance)) {
		std::fprintf(stderr, "%s: got %.20Lg expected %.20Lg error %.4Lg\n",
				name, got, expected, error);
		std::exit(1);
	}
}

int main()
{
	/* Rectangle [-1,1] x [-2,2], centered at the generator. */
	const double rx[] = {-1.0, 1.0, 1.0, -1.0};
	const double ry[] = {-2.0, -2.0, 2.0, 2.0};
	GfsPolygonMomentsReference rectangle{};
	if (!gfs_polygon_moments_reference(rx, ry, 4, &rectangle)) return 1;
	require_close("rectangle area", rectangle.area, 8.0L, 1e-18L);
	require_close("rectangle cx", rectangle.centroid_x, 0.0L, 1e-18L);
	require_close("rectangle cy", rectangle.centroid_y, 0.0L, 1e-18L);
	require_close("rectangle Qxx", rectangle.qxx_area, 8.0L / 3.0L, 1e-18L);
	require_close("rectangle Qxy", rectangle.qxy_area, 0.0L, 1e-18L);
	require_close("rectangle Qyy", rectangle.qyy_area, 32.0L / 3.0L, 1e-18L);

	/* Right triangle (0,0),(2,0),(0,1), with a nonzero product moment. */
	const double tx[] = {0.0, 2.0, 0.0};
	const double ty[] = {0.0, 0.0, 1.0};
	GfsPolygonMomentsReference triangle{};
	if (!gfs_polygon_moments_reference(tx, ty, 3, &triangle)) return 1;
	require_close("triangle area", triangle.area, 1.0L, 1e-18L);
	require_close("triangle cx", triangle.centroid_x, 2.0L / 3.0L, 1e-18L);
	require_close("triangle cy", triangle.centroid_y, 1.0L / 3.0L, 1e-18L);
	require_close("triangle Qxx", triangle.qxx_area, 2.0L / 9.0L, 1e-18L);
	require_close("triangle Qxy", triangle.qxy_area, -1.0L / 18.0L, 1e-18L);
	require_close("triangle Qyy", triangle.qyy_area, 1.0L / 18.0L, 1e-18L);

	/* A local translation must shift only the centroid, not central Q. */
	const double tx_shift[] = {0.375, 2.375, 0.375};
	const double ty_shift[] = {-0.25, -0.25, 0.75};
	GfsPolygonMomentsReference shifted{};
	if (!gfs_polygon_moments_reference(tx_shift, ty_shift, 3, &shifted)) return 1;
	require_close("shifted cx", shifted.centroid_x,
			triangle.centroid_x + 0.375L, 1e-18L);
	require_close("shifted cy", shifted.centroid_y,
			triangle.centroid_y - 0.25L, 1e-18L);
	require_close("shifted Qxx", shifted.qxx_area, triangle.qxx_area, 1e-18L);
	require_close("shifted Qxy", shifted.qxy_area, triangle.qxy_area, 1e-18L);
	require_close("shifted Qyy", shifted.qyy_area, triangle.qyy_area, 1e-18L);

	/* The production clipping convention is CCW; reject reversed and flat data. */
	const double cw_x[] = {0.0, 0.0, 2.0};
	const double cw_y[] = {0.0, 1.0, 0.0};
	GfsPolygonMomentsReference invalid{};
	if (gfs_polygon_moments_reference(cw_x, cw_y, 3, &invalid)) return 1;
	const double flat_x[] = {0.0, 1.0, 2.0};
	const double flat_y[] = {0.0, 0.0, 0.0};
	if (gfs_polygon_moments_reference(flat_x, flat_y, 3, &invalid)) return 1;

	std::puts("PASS: signed polygon area, centroid and central Q reference");
	return 0;
}
