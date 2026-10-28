#ifndef GFS_SPIN_FACE_EXCHANGE_20261027_H
#define GFS_SPIN_FACE_EXCHANGE_20261027_H

#include <math.h>

#if defined(__CUDACC__)
#define GFS_SPIN_HD __host__ __device__
#else
#define GFS_SPIN_HD
#endif

/* Algebraic two-cell face-exchange ledger. x,y are mass centroids, not mesh
 * generators; spin and inertia are about those centroids. Energy includes
 * bulk, rotational and thermal terms. This isolated flux substep transports
 * the first mass moment through the face centroid. The caller supplies
 * physical fluxes and post-exchange inertias. It does not compute a Riemann
 * flux, polygon moments, centroid volume terms, generator motion, or RK4. */
struct gfs_spin_cell {
    double x,y,m,px,py,energy,spin,inertia;
};

struct gfs_spin_face_transfer {
    double cx,cy,dm,dpx,dpy,denergy,inertia_i_new,inertia_j_new;
};

GFS_SPIN_HD static inline bool gfs_spin_cell_finite(const gfs_spin_cell &a) {
    return isfinite(a.x) && isfinite(a.y) && isfinite(a.m) &&
           isfinite(a.px) && isfinite(a.py) && isfinite(a.energy) &&
           isfinite(a.spin) && isfinite(a.inertia);
}

GFS_SPIN_HD static inline double gfs_spin_thermal_energy(const gfs_spin_cell &a) {
    return a.energy-(a.px*a.px+a.py*a.py)/(2.*a.m)
                   -a.spin*a.spin/(2.*a.inertia);
}

/* Positive dm, dp and denergy flow from i to j. Return 0 only after both
 * updated cells are finite and have positive mass, inertia and thermal
 * energy. Failure leaves both input cells bitwise unchanged. */
GFS_SPIN_HD static inline int gfs_spin_face_exchange(
    gfs_spin_cell *i, gfs_spin_cell *j, const gfs_spin_face_transfer &f)
{
    if (!i || !j || i==j || !gfs_spin_cell_finite(*i) ||
        !gfs_spin_cell_finite(*j) || !isfinite(f.cx) || !isfinite(f.cy) ||
        !isfinite(f.dm) || !isfinite(f.dpx) || !isfinite(f.dpy) ||
        !isfinite(f.denergy) || !isfinite(f.inertia_i_new) ||
        !isfinite(f.inertia_j_new)) return -1;
    if (!(i->m>0.) || !(j->m>0.) || !(i->inertia>0.) || !(j->inertia>0.) ||
        !(f.inertia_i_new>0.) || !(f.inertia_j_new>0.)) return -2;
    if (!(gfs_spin_thermal_energy(*i)>0.) || !(gfs_spin_thermal_energy(*j)>0.)) return -3;
    gfs_spin_cell a=*i,b=*j;
    a.m-=f.dm; b.m+=f.dm;
    if (!isfinite(a.m) || !isfinite(b.m)) return -1;
    if (!(a.m>0.) || !(b.m>0.)) return -2;
    a.px-=f.dpx; a.py-=f.dpy;
    b.px+=f.dpx; b.py+=f.dpy;
    a.energy-=f.denergy; b.energy+=f.denergy;
    a.x+=f.dm*(i->x-f.cx)/a.m; a.y+=f.dm*(i->y-f.cy)/a.m;
    b.x+=f.dm*(f.cx-j->x)/b.m; b.y+=f.dm*(f.cy-j->y)/b.m;
    a.spin-=((f.cx-i->x)*f.dpy-(f.cy-i->y)*f.dpx);
    b.spin+=((f.cx-j->x)*f.dpy-(f.cy-j->y)*f.dpx);
    a.spin+=(i->x-a.x)*a.py-(i->y-a.y)*a.px;
    b.spin+=(j->x-b.x)*b.py-(j->y-b.y)*b.px;
    a.inertia=f.inertia_i_new; b.inertia=f.inertia_j_new;
    if (!gfs_spin_cell_finite(a) || !gfs_spin_cell_finite(b)) return -1;
    if (!(gfs_spin_thermal_energy(a)>0.) || !(gfs_spin_thermal_energy(b)>0.)) return -3;
    *i=a; *j=b;
    return 0;
}

#undef GFS_SPIN_HD
#endif
