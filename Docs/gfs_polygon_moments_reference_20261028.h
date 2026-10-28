#ifndef GFS_POLYGON_MOMENTS_REFERENCE_20261028_H
#define GFS_POLYGON_MOMENTS_REFERENCE_20261028_H

#include <cmath>
#include <cstddef>

struct GfsPolygonMomentsReference {
	long double area;
	long double centroid_x;
	long double centroid_y;
	long double qxx_area;
	long double qxy_area;
	long double qyy_area;
};

/* Independent signed polygon sums for CCW vertices in generator-relative
 * coordinates. Q components are area moments about the polygon centroid. */
inline bool gfs_polygon_moments_reference(
		const double *x, const double *y, std::size_t n,
		GfsPolygonMomentsReference *out)
{
	if (!x || !y || !out || n < 3) return false;
	long double area2 = 0.0L;
	long double cx_sum = 0.0L, cy_sum = 0.0L;
	long double jxx12 = 0.0L, jyy12 = 0.0L, jxy24 = 0.0L;
	for (std::size_t i = 0; i < n; ++i) {
		const std::size_t j = (i + 1) % n;
		const long double xi = x[i], yi = y[i];
		const long double xj = x[j], yj = y[j];
		const long double cross = xi * yj - xj * yi;
		area2 += cross;
		cx_sum += (xi + xj) * cross;
		cy_sum += (yi + yj) * cross;
		jxx12 += (xi * xi + xi * xj + xj * xj) * cross;
		jyy12 += (yi * yi + yi * yj + yj * yj) * cross;
		jxy24 += (2.0L * xi * yi + xi * yj + xj * yi + 2.0L * xj * yj) * cross;
	}
	if (!(area2 > 0.0L) || !std::isfinite(area2)) return false;
	const long double area = 0.5L * area2;
	const long double cx = cx_sum / (3.0L * area2);
	const long double cy = cy_sum / (3.0L * area2);
	const long double jxx = jxx12 / 12.0L;
	const long double jyy = jyy12 / 12.0L;
	const long double jxy = jxy24 / 24.0L;
	const long double qxx = jxx - area * cx * cx;
	const long double qyy = jyy - area * cy * cy;
	const long double qxy = jxy - area * cx * cy;
	const long double det = qxx * qyy - qxy * qxy;
	if (!std::isfinite(cx) || !std::isfinite(cy) ||
			!std::isfinite(qxx) || !std::isfinite(qxy) ||
			!std::isfinite(qyy) || !(qxx > 0.0L) || !(qyy > 0.0L) ||
			!(det > 0.0L)) return false;
	out->area = area;
	out->centroid_x = cx;
	out->centroid_y = cy;
	out->qxx_area = qxx;
	out->qxy_area = qxy;
	out->qyy_area = qyy;
	return true;
}

#endif
