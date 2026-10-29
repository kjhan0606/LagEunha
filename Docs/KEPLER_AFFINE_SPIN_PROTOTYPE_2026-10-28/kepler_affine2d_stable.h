#ifndef KEPLER_AFFINE2D_H
#define KEPLER_AFFINE2D_H

/* Exact 2-D polygon moments and conservative affine-velocity recovery.
 * Quantities are extensive and use double precision independently of the
 * solver's PosType. This header has no MPI/CUDA dependencies. */

#include <math.h>
#include <stddef.h>

typedef struct KeplerPolygonMoments2D {
	double area;
	double cx, cy;
	double sxx, sxy, syy; /* central area moments */
} KeplerPolygonMoments2D;

typedef struct KeplerAffineConserved2D {
	double mass;
	double hx, hy;                 /* integral rho x_a */
	double px, py;                 /* integral rho u_b */
	double qxx, qxy, qyy;         /* integral rho x_a x_b */
	double rxx, rxy, ryx, ryy;    /* integral rho x_a u_b */
	double etotal;                /* thermal + all resolved kinetic energy */
} KeplerAffineConserved2D;

typedef struct KeplerAffineDerived2D {
	double cx, cy;
	double sxx, sxy, syy;         /* central mass moments */
	double gxx, gxy, gyx, gyy;    /* u = P/M + G (x-C) */
	double angular_momentum_z;
	double ebulk, eaffine, ethermal;
} KeplerAffineDerived2D;

/* Lagrangian fixed-mass state about its moving mass centroid. S is the
 * central mass covariance and B_ab = integral rho r_a (u_b-U_b) dA. */
typedef struct KeplerAffineStateLagrangian2D {
	double cx, cy;
	double sxx, sxy, syy;
	double bxx, bxy, byx, byy;
	double etotal;
} KeplerAffineStateLagrangian2D;

typedef struct KeplerAffineRateLagrangian2D {
	double dcx, dcy;
	double dsxx, dsxy, dsyy;
	double dbxx, dbxy, dbyx, dbyy;
	double dpx, dpy; /* extensive momentum rate */
	double detotal;
} KeplerAffineRateLagrangian2D;

typedef struct KeplerAffineFace2D {
	double length;
	double nx, ny;       /* outward normal from left cell to right cell */
	double rlx, rly;     /* face centroid relative to each mass centroid */
	double rrx, rry;
	double ux, uy;       /* one shared face velocity */
	double pressure;     /* one shared normal traction */
} KeplerAffineFace2D;

/* Polygon integrals are exact for straight edges. Either winding direction
 * is accepted; returned area and central moments are positive-orientation. */
static inline int kepler_polygon_moments2d(const double *x, const double *y,
		size_t n, KeplerPolygonMoments2D *out)
{
	double twice_area = 0.0, cx6 = 0.0, cy6 = 0.0;
	double raw_xx12 = 0.0, raw_xy24 = 0.0, raw_yy12 = 0.0;
	if (!x || !y || !out || n < 3) return 0;
	/* Compute central moments in a local frame. Global raw moments minus
	 * area*centroid^2 lose many digits for small, translated cells. */
	const double xref = x[0], yref = y[0];
	if (!isfinite(xref) || !isfinite(yref)) return 0;
	for (size_t i = 0; i < n; ++i) {
		size_t j = (i + 1) % n;
		double xi = x[i]-xref, yi = y[i]-yref;
		double xj = x[j]-xref, yj = y[j]-yref;
		double cross = xi*yj - xj*yi;
		if (!isfinite(cross)) return 0;
		twice_area += cross;
		cx6 += (xi + xj)*cross;
		cy6 += (yi + yj)*cross;
		raw_xx12 += (xi*xi + xi*xj + xj*xj)*cross;
		raw_yy12 += (yi*yi + yi*yj + yj*yj)*cross;
		raw_xy24 += (2.0*xi*yi + xi*yj + xj*yi + 2.0*xj*yj)*cross;
	}
	if (!isfinite(twice_area) || fabs(twice_area) <= 1.0e-14) return 0;
	double orient = twice_area > 0.0 ? 1.0 : -1.0;
	double area = 0.5*fabs(twice_area);
	double cx = orient*cx6/(3.0*fabs(twice_area));
	double cy = orient*cy6/(3.0*fabs(twice_area));
	double raw_xx = orient*raw_xx12/12.0;
	double raw_xy = orient*raw_xy24/24.0;
	double raw_yy = orient*raw_yy12/12.0;
	out->area = area;
	out->cx = xref + cx;
	out->cy = yref + cy;
	out->sxx = raw_xx - area*cx*cx;
	out->sxy = raw_xy - area*cx*cy;
	out->syy = raw_yy - area*cy*cy;
	return isfinite(out->sxx) && isfinite(out->sxy) && isfinite(out->syy);
}

/* Exact first and second line moments on a straight face. */
static inline int kepler_segment_moments2d(double ax, double ay,
		double bx, double by, double *length,
		double first[2], double second[3])
{
	double dx = bx - ax, dy = by - ay;
	double len = hypot(dx, dy);
	double mx = 0.5*(ax + bx), my = 0.5*(ay + by);
	if (!(len > 0.0) || !isfinite(len) || !first || !second) return 0;
	if (length) *length = len;
	first[0] = len*mx;
	first[1] = len*my;
	second[0] = len*(mx*mx + dx*dx/12.0);
	second[1] = len*(mx*my + dx*dy/12.0);
	second[2] = len*(my*my + dy*dy/12.0);
	return 1;
}

/* Recover the affine velocity field and energy split from extensive moments.
 * Returns 1 on a valid positive-definite covariance and finite state. */
static inline int kepler_affine_derive2d(const KeplerAffineConserved2D *u,
		KeplerAffineDerived2D *d)
{
	if (!u || !d || !(u->mass > 0.0) || !isfinite(u->mass)) return 0;
	d->cx = u->hx/u->mass;
	d->cy = u->hy/u->mass;
	d->sxx = u->qxx - u->hx*u->hx/u->mass;
	d->sxy = u->qxy - u->hx*u->hy/u->mass;
	d->syy = u->qyy - u->hy*u->hy/u->mass;
	double det = d->sxx*d->syy - d->sxy*d->sxy;
	if (!(d->sxx > 0.0) || !(d->syy > 0.0) || !(det > 0.0) || !isfinite(det)) return 0;
	const double ixx = d->syy/det, ixy = -d->sxy/det, iyy = d->sxx/det;
	double bxx = u->rxx - u->hx*u->px/u->mass;
	double bxy = u->rxy - u->hx*u->py/u->mass;
	double byx = u->ryx - u->hy*u->px/u->mass;
	double byy = u->ryy - u->hy*u->py/u->mass;
	/* G = B^T S^{-1}; B_{a b}=integral rho (x_a-C_a)(u_b-U_b). */
	d->gxx = bxx*ixx + byx*ixy;
	d->gxy = bxx*ixy + byx*iyy;
	d->gyx = bxy*ixx + byy*ixy;
	d->gyy = bxy*ixy + byy*iyy;
	d->angular_momentum_z = bxy - byx;
	d->ebulk = 0.5*(u->px*u->px + u->py*u->py)/u->mass;
	d->eaffine = 0.5*(
		d->gxx*d->gxx*d->sxx + 2.0*d->gxx*d->gxy*d->sxy + d->gxy*d->gxy*d->syy +
		d->gyx*d->gyx*d->sxx + 2.0*d->gyx*d->gyy*d->sxy + d->gyy*d->gyy*d->syy);
	d->ethermal = u->etotal - d->ebulk - d->eaffine;
	return isfinite(d->cx) && isfinite(d->cy) && isfinite(d->gxx) &&
		isfinite(d->gxy) && isfinite(d->gyx) && isfinite(d->gyy) &&
		isfinite(d->angular_momentum_z) && isfinite(d->ethermal);
}

static inline int kepler_affine_state_from_polygon2d(
		const KeplerPolygonMoments2D *p, double mass,
		double vx, double vy, double gxx, double gxy,
		double gyx, double gyy, double ethermal,
		KeplerAffineStateLagrangian2D *s)
{
	if (!p || !s || !(p->area > 0.0) || !(mass > 0.0) || ethermal < 0.0) return 0;
	double rho = mass/p->area;
	s->cx = p->cx; s->cy = p->cy;
	s->sxx = rho*p->sxx; s->sxy = rho*p->sxy; s->syy = rho*p->syy;
	/* B = S G^T. */
	s->bxx = s->sxx*gxx + s->sxy*gxy;
	s->bxy = s->sxx*gyx + s->sxy*gyy;
	s->byx = s->sxy*gxx + s->syy*gxy;
	s->byy = s->sxy*gyx + s->syy*gyy;
	double eaff = 0.5*(
		gxx*gxx*s->sxx + 2.0*gxx*gxy*s->sxy + gxy*gxy*s->syy +
		gyx*gyx*s->sxx + 2.0*gyx*gyy*s->sxy + gyy*gyy*s->syy);
	s->etotal = ethermal + 0.5*mass*(vx*vx + vy*vy) + eaff;
	return isfinite(s->etotal);
}

/* Exact gradient of the Kepler IC velocity used by kp_mkinitial:
 * u=(-dy,dx) (r^2+eps^2)^(-3/4). Coordinates are relative to the central
 * mass; eps is the softening length (zero gives the point-mass disk). */
static inline int kepler_softened_velocity_gradient2d(
		double dx, double dy, double eps,
		double *gxx, double *gxy, double *gyx, double *gyy)
{
	double r2 = dx*dx + dy*dy;
	double s = r2 + eps*eps;
	if (!(s > 0.0) || !isfinite(s) || !gxx || !gxy || !gyx || !gyy) return 0;
	double f = pow(s, -0.75);
	double q = -1.5*pow(s, -1.75); /* (1/r) df/dr */
	*gxx = -dx*dy*q;
	*gxy = -f - dy*dy*q;
	*gyx =  f + dx*dx*q;
	*gyy =  dx*dy*q;
	return isfinite(*gxx) && isfinite(*gxy) && isfinite(*gyx) && isfinite(*gyy);
}

/* Jacobian of the softened unit-mass Kepler acceleration
 * a=-r/(r^2+eps^2)^(3/2). The row index is acceleration component and the
 * column index is spatial derivative: A_ij = d a_i / d x_j. This tidal
 * tensor enters the affine cell-moment source and is evaluated at the cell's
 * mass centroid, not at the mesh generator. */
static inline int kepler_softened_accel_gradient2d(
		double dx, double dy, double eps,
		double *axx, double *axy, double *ayx, double *ayy)
{
	double s = dx*dx + dy*dy + eps*eps;
	if (!(s > 0.0) || !isfinite(s) || !axx || !axy || !ayx || !ayy) return 0;
	double s_m3 = pow(s, -1.5);
	double s_m5 = pow(s, -2.5);
	*axx = -s_m3 + 3.0*dx*dx*s_m5;
	*axy =  3.0*dx*dy*s_m5;
	*ayx =  *axy;
	*ayy = -s_m3 + 3.0*dy*dy*s_m5;
	return isfinite(*axx) && isfinite(*axy) && isfinite(*ayx) && isfinite(*ayy);
}

static inline int kepler_affine_state_derive2d(
		const KeplerAffineStateLagrangian2D *s, double mass,
		double vx, double vy, KeplerAffineDerived2D *d)
{
	if (!s || !d || !(mass > 0.0)) return 0;
	d->cx = s->cx; d->cy = s->cy;
	d->sxx = s->sxx; d->sxy = s->sxy; d->syy = s->syy;
	double det = s->sxx*s->syy - s->sxy*s->sxy;
	if (!(s->sxx > 0.0) || !(s->syy > 0.0) || !(det > 0.0) || !isfinite(det)) return 0;
	double ixx = s->syy/det, ixy = -s->sxy/det, iyy = s->sxx/det;
	d->gxx = s->bxx*ixx + s->byx*ixy;
	d->gxy = s->bxx*ixy + s->byx*iyy;
	d->gyx = s->bxy*ixx + s->byy*ixy;
	d->gyy = s->bxy*ixy + s->byy*iyy;
	d->angular_momentum_z = s->bxy - s->byx;
	d->ebulk = 0.5*mass*(vx*vx + vy*vy);
	d->eaffine = 0.5*(
		d->gxx*d->gxx*s->sxx + 2.0*d->gxx*d->gxy*s->sxy + d->gxy*d->gxy*s->syy +
		d->gyx*d->gyx*s->sxx + 2.0*d->gyx*d->gyy*s->sxy + d->gyy*d->gyy*s->syy);
	d->ethermal = s->etotal - d->ebulk - d->eaffine;
	return isfinite(d->gxx) && isfinite(d->gxy) && isfinite(d->gyx) &&
		isfinite(d->gyy) && isfinite(d->angular_momentum_z) &&
		isfinite(d->ethermal);
}

/* Add one pressure face exactly once to both cells. The routine implements
 * action/reaction force, the matching pressure work, and intrinsic moment
 * flux. Body forces and the volume stress source are accumulated separately. */
static inline int kepler_affine_add_shared_face2d(
		KeplerAffineRateLagrangian2D *left,
		KeplerAffineRateLagrangian2D *right,
		const KeplerAffineFace2D *f)
{
	if (!left || !f || !(f->length > 0.0) || f->pressure < 0.0) return 0;
	double fx = -f->pressure*f->length*f->nx;
	double fy = -f->pressure*f->length*f->ny;
	left->dpx += fx; left->dpy += fy;
	left->dbxx += f->rlx*fx;
	left->dbxy += f->rlx*fy;
	left->dbyx += f->rly*fx;
	left->dbyy += f->rly*fy;
	left->detotal -= f->pressure*f->length*(f->ux*f->nx + f->uy*f->ny);
	if (right) {
		right->dpx -= fx; right->dpy -= fy;
		right->dbxx -= f->rrx*fx;
		right->dbxy -= f->rrx*fy;
		right->dbyx -= f->rry*fx;
		right->dbyy -= f->rry*fy;
		right->detotal += f->pressure*f->length*(f->ux*f->nx + f->uy*f->ny);
	}
	return 1;
}

/* Fixed-mass cell volume and body-force sources. Agrad is the spatial
 * acceleration gradient, with rows as acceleration components. */
static inline int kepler_affine_add_cell_sources2d(
		KeplerAffineRateLagrangian2D *r,
		const KeplerAffineStateLagrangian2D *s,
		double mass, double vx, double vy, double gamma,
		double ax, double ay, double axx, double axy,
		double ayx, double ayy)
{
	KeplerAffineDerived2D d;
	if (!r || !kepler_affine_state_derive2d(s, mass, vx, vy, &d) ||
	    !(gamma > 1.0) || !(d.ethermal >= 0.0)) return 0;
	r->dcx += vx; r->dcy += vy;
	r->dsxx += 2.0*s->bxx;
	r->dsxy += s->bxy + s->byx;
	r->dsyy += 2.0*s->byy;
	double wxx = d.gxx*d.gxx*s->sxx + 2.0*d.gxx*d.gxy*s->sxy + d.gxy*d.gxy*s->syy;
	double wxy = d.gxx*d.gyx*s->sxx + (d.gxx*d.gyy + d.gxy*d.gyx)*s->sxy + d.gxy*d.gyy*s->syy;
	double wyy = d.gyx*d.gyx*s->sxx + 2.0*d.gyx*d.gyy*s->sxy + d.gyy*d.gyy*s->syy;
	double pvolume = (gamma - 1.0)*d.ethermal;
	r->dbxx += wxx + pvolume + s->sxx*axx + s->sxy*axy;
	r->dbxy += wxy + s->sxx*ayx + s->sxy*ayy;
	r->dbyx += wxy + s->sxy*axx + s->syy*axy;
	r->dbyy += wyy + pvolume + s->sxy*ayx + s->syy*ayy;
	r->dpx += mass*ax; r->dpy += mass*ay;
	r->detotal += mass*(vx*ax + vy*ay) +
		s->bxx*axx + s->bxy*ayx + s->byx*axy + s->byy*ayy;
	return 1;
}

#endif /* KEPLER_AFFINE2D_H */
