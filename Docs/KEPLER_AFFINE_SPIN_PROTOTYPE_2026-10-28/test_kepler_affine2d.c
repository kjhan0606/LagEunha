#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "kepler_affine2d.h"

static int close_to(double a, double b) { return fabs(a-b) < 1.0e-12; }

int main(void)
{
	double x[4] = {-1, 1, 1, -1}, y[4] = {-1, -1, 1, 1};
	KeplerPolygonMoments2D p;
	assert(kepler_polygon_moments2d(x, y, 4, &p));
	assert(close_to(p.area, 4.0));
	assert(close_to(p.cx, 0.0) && close_to(p.cy, 0.0));
	assert(close_to(p.sxx, 4.0/3.0) && close_to(p.sxy, 0.0));
	assert(close_to(p.syy, 4.0/3.0));
	for (int i = 0; i < 2; ++i) {
		double t = x[i]; x[i] = x[3-i]; x[3-i] = t;
		t = y[i]; y[i] = y[3-i]; y[3-i] = t;
	}
	assert(kepler_polygon_moments2d(x, y, 4, &p));
	assert(close_to(p.area, 4.0) && close_to(p.sxx, 4.0/3.0));

	KeplerAffineConserved2D u = {0};
	KeplerAffineDerived2D d;
	u.mass = 4.0; u.qxx = u.qyy = 4.0/3.0;
	/* Uniform square, rigid rotation Omega=3, plus thermal energy 5. */
	u.rxy = 4.0; u.ryx = -4.0; u.etotal = 17.0;
	assert(kepler_affine_derive2d(&u, &d));
	assert(close_to(d.gxx, 0.0) && close_to(d.gxy, -3.0));
	assert(close_to(d.gyx, 3.0) && close_to(d.gyy, 0.0));
	assert(close_to(d.angular_momentum_z, 8.0));
	assert(close_to(d.eaffine, 12.0) && close_to(d.ethermal, 5.0));

	KeplerAffineStateLagrangian2D s;
	assert(kepler_affine_state_from_polygon2d(&p, 4.0, 0, 0,
		0, -3, 3, 0, 5.0, &s));
	assert(kepler_affine_state_derive2d(&s, 4.0, 0, 0, &d));
	assert(close_to(d.gxy, -3.0) && close_to(d.gyx, 3.0));
	assert(close_to(d.angular_momentum_z, 8.0));
	assert(close_to(d.eaffine, 12.0) && close_to(d.ethermal, 5.0));
	KeplerAffineRateLagrangian2D cell_rate = {0};
	assert(kepler_affine_add_cell_sources2d(&cell_rate, &s, 4.0, 0, 0,
		1.4, 0, 0, 0, 0, 0, 0));
	assert(close_to(cell_rate.dsxx, 0.0) && close_to(cell_rate.dsxy, 0.0));
	assert(close_to(cell_rate.dsyy, 0.0));
	assert(close_to(cell_rate.dbxy, cell_rate.dbyx));
	assert(close_to(cell_rate.detotal, 0.0));

	KeplerAffineRateLagrangian2D left = {0}, right = {0};
	KeplerAffineFace2D face = {
		1.0, 1.0, 0.0, 1.0, 1.0, -1.0, 1.0, 0.5, 0.2, 2.0
	};
	assert(kepler_affine_add_shared_face2d(&left, &right, &face));
	assert(close_to(left.dpx + right.dpx, 0.0));
	assert(close_to(left.dpy + right.dpy, 0.0));
	assert(close_to(left.detotal + right.detotal, 0.0));
	assert(close_to(left.dbxy-left.dbyx, 2.0));
	assert(close_to(right.dbxy-right.dbyx, -2.0));

	double gxx, gxy, gyx, gyy;
	double dx = 0.8, dy = -0.4, eps = 0.01175, h = 1.0e-6;
	assert(kepler_softened_velocity_gradient2d(dx, dy, eps,
		&gxx, &gxy, &gyx, &gyy));
	double f0 = pow(dx*dx + dy*dy + eps*eps, -0.75);
	double fpx = pow((dx+h)*(dx+h) + dy*dy + eps*eps, -0.75);
	double fmx = pow((dx-h)*(dx-h) + dy*dy + eps*eps, -0.75);
	double fpy = pow(dx*dx + (dy+h)*(dy+h) + eps*eps, -0.75);
	double fmy = pow(dx*dx + (dy-h)*(dy-h) + eps*eps, -0.75);
	assert(fabs(gxx - (-dy*fpx + dy*fmx)/(2*h)) < 1.0e-8);
	assert(fabs(gyx - ((dx+h)*fpx - (dx-h)*fmx)/(2*h)) < 1.0e-8);
	assert(fabs(gxy - (-(dy+h)*fpy + (dy-h)*fmy)/(2*h)) < 1.0e-8);
	assert(fabs(gyy - (dx*fpy - dx*fmy)/(2*h)) < 1.0e-8);
	assert(f0 > 0.0);
	double axx, axy, ayx, ayy;
	assert(kepler_softened_accel_gradient2d(dx, dy, eps,
		&axx, &axy, &ayx, &ayy));
	double apx = -(dx+h)*pow((dx+h)*(dx+h)+dy*dy+eps*eps, -1.5);
	double amx = -(dx-h)*pow((dx-h)*(dx-h)+dy*dy+eps*eps, -1.5);
	double axyp = -dx*pow(dx*dx+(dy+h)*(dy+h)+eps*eps, -1.5);
	double axym = -dx*pow(dx*dx+(dy-h)*(dy-h)+eps*eps, -1.5);
	double ayxp = -dy*pow((dx+h)*(dx+h)+dy*dy+eps*eps, -1.5);
	double ayxm = -dy*pow((dx-h)*(dx-h)+dy*dy+eps*eps, -1.5);
	double ayyp = -(dy+h)*pow(dx*dx+(dy+h)*(dy+h)+eps*eps, -1.5);
	double ayym = -(dy-h)*pow(dx*dx+(dy-h)*(dy-h)+eps*eps, -1.5);
	assert(fabs(axx - (apx-amx)/(2*h)) < 1.0e-8);
	assert(fabs(axy - (axyp-axym)/(2*h)) < 1.0e-8);
	assert(fabs(ayx - (ayxp-ayxm)/(2*h)) < 1.0e-8);
	assert(fabs(ayy - (ayyp-ayym)/(2*h)) < 1.0e-8);

	double length, first[2], second[3];
	assert(kepler_segment_moments2d(0, 0, 3, 0, &length, first, second));
	assert(close_to(length, 3.0) && close_to(first[0], 4.5));
	assert(close_to(second[0], 9.0) && close_to(second[1], 0.0));
	puts("kepler_affine2d analytic checks: PASS");
	return 0;
}
